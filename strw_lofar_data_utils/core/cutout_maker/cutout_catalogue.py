from pandas.core.frame import DataFrame
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

from .source_blob import SourceBlob
from .cutout import Cutout


class CutoutCatalogue:
    _catalog_cache: dict[tuple[int, str, str, str, str], dict] = {}

    def __init__(
            self,
            catalogue: DataFrame,
            cutout: Cutout,
            ra_col: str = "RA",
            dec_col: str = "DEC",
            comp_col: str = "Component_Name",
            source_col: str = "Parent_Source",
            ):
        """
        Initialize a CutoutCatalogue object that constrains a catalogue to objects within a cutout area.
        The catalogue is expected to have columns for RA, Dec, Component Name, and Source Name, otherwise
        a ValueError is raised.
        
        :param catalogue: Pandas DataFrame containing the source catalogue
        :type catalogue: DataFrame
        :param cutout: Cutout object defining the cutout area
        :type cutout: Cutout
        :param ra_col: Name of the RA column in the catalogue (default: "RA")
        :type ra_col: str
        :param dec_col: Name of the Dec column in the catalogue (default: "DEC")
        :type dec_col: str
        :param comp_col: Name of the Component Name column in the catalogue (default: "Component_Name").
        Optional - if your catalogue does not have a component column, set comp_col=None when initializing CutoutCatalogue.
        :type comp_col: str
        :param source_col: Name of the Source Name column in the catalogue (default: "Parent_Source"). This
        column should be set to the name of the column in the PyBDSF-style catalogue that contains the
        source names (which can have multiple components). This is used to group components into sources.
        If your catalogue does not have a source column each component is treated as a separate source.
        :type source_col: str
        """
        # Make sure catalogue has required columns
        for col in [ra_col, dec_col, source_col]:
            if col not in catalogue.columns:
                raise ValueError(f"Catalogue is missing required column: {col}")
            
        if comp_col is not None:
            if comp_col not in catalogue.columns:
                raise ValueError(
                    f"Catalogue is missing component column: {comp_col}. "
                    "If your catalogue does not have a component column, set comp_col=None when initializing CutoutCatalogue."
                )
            
        # Store cutout properties
        self.cutout_shape = (cutout.size_pixels, cutout.size_pixels)
        self.cutout_wcs = cutout.get_wcs()
        self.source_col = source_col
        self.ra_col = ra_col
        self.dec_col = dec_col
        self.comp_col = comp_col

        cache_key = (id(catalogue), ra_col, dec_col, comp_col, source_col)
        if cache_key not in self._catalog_cache:
            ra_all = np.asarray(catalogue[ra_col].to_numpy(), dtype=np.float64)
            dec_all = np.asarray(catalogue[dec_col].to_numpy(), dtype=np.float64)
            sort_idx = np.argsort(ra_all)
            self._catalog_cache[cache_key] = {
                "ra_sorted": ra_all[sort_idx],
                "dec_sorted": dec_all[sort_idx],
                "source_sorted": catalogue[source_col].to_numpy()[sort_idx],
                "comp_sorted": catalogue[comp_col].to_numpy()[sort_idx] if comp_col is not None else None,
                "sort_idx": sort_idx,
            }
        self._cache = self._catalog_cache[cache_key]
        
        # Constrain catalogue to objects within cutout area
        size_arcmin = cutout.size_arcmin
        cutout_ra = cutout.ra
        cutout_dec = cutout.dec
        self.constrained_catalogue = self._constrain_catalogue(
            catalogue,
            size_arcmin,
            cutout_ra,
            cutout_dec,
            ra_col,
            dec_col
        )

        # Precompute cutout inclusion once per instance and reuse in all retrieval methods
        if len(self.constrained_catalogue) == 0:
            self._filtered_ra = np.array([], dtype=np.float64)
            self._filtered_dec = np.array([], dtype=np.float64)
            self._filtered_source = np.array([], dtype=object)
            self._filtered_comp = np.array([], dtype=object) if comp_col is not None else None
            self._filtered_x = np.array([], dtype=np.int32)
            self._filtered_y = np.array([], dtype=np.int32)
        else:
            constrained_ra = np.asarray(self.constrained_catalogue[self.ra_col].to_numpy(), dtype=np.float64)
            constrained_dec = np.asarray(self.constrained_catalogue[self.dec_col].to_numpy(), dtype=np.float64)
            constrained_source = self.constrained_catalogue[self.source_col].to_numpy()
            constrained_comp = self.constrained_catalogue[self.comp_col].to_numpy() if self.comp_col is not None else None

            coords = SkyCoord(ra=constrained_ra * u.deg, dec=constrained_dec * u.deg, frame='icrs')
            x_pixel, y_pixel = self.cutout_wcs.world_to_pixel(coords)
            x_pixel = np.rint(x_pixel).astype(np.int32)
            y_pixel = np.rint(y_pixel).astype(np.int32)

            in_cutout_mask = (
                (x_pixel >= 0) & (x_pixel < self.cutout_shape[1]) &
                (y_pixel >= 0) & (y_pixel < self.cutout_shape[0])
            )

            self._filtered_ra = constrained_ra[in_cutout_mask]
            self._filtered_dec = constrained_dec[in_cutout_mask]
            self._filtered_source = constrained_source[in_cutout_mask]
            self._filtered_comp = constrained_comp[in_cutout_mask] if constrained_comp is not None else None
            self._filtered_x = x_pixel[in_cutout_mask]
            self._filtered_y = y_pixel[in_cutout_mask]

    def get_source_blobs_from_catalogue(self, unique_objects: bool = True) -> dict[str, SourceBlob]:
        """
        Create SourceBlob instances for each unique object in the constrained catalogue
        that fall within the cutout area. Returns a dictionary mapping source names
        to SourceBlob instances. This method checks if the objects are within the cutout
        area using the cutout WCS and shape and the check_if_in_cutout method of SourceBlob.
        """
        # Get unique source names
        if unique_objects:
            # We use only the sources (which can have multiple components)
            target_column = self.source_col
        else:
            # We use the components (which can be multiple per source)
            target_column = self.comp_col

        # Keep ordering consistent with original implementation (.unique() order on constrained catalogue)
        unique_sources = self.constrained_catalogue[target_column].unique()

        if len(unique_sources) == 0 or len(self._filtered_ra) == 0:
            return {}

        filtered_target_values = self._filtered_source if unique_objects else self._filtered_comp

        # Aggregate filtered rows by target object in a single pass
        grouped = {}
        for key, ra, dec, x_pos, y_pos in zip(
            filtered_target_values,
            self._filtered_ra,
            self._filtered_dec,
            self._filtered_x,
            self._filtered_y,
        ):
            if key not in grouped:
                grouped[key] = {
                    "ra": [],
                    "dec": [],
                    "x": [],
                    "y": [],
                }
            grouped[key]["ra"].append(float(ra))
            grouped[key]["dec"].append(float(dec))
            grouped[key]["x"].append(int(x_pos))
            grouped[key]["y"].append(int(y_pos))

        source_blobs = {}
        for source in unique_sources:
            if source not in grouped:
                continue

            group = grouped[source]
            sb = SourceBlob(group["ra"], group["dec"])
            sb.set_precomputed_positions(group["x"], group["y"])
            source_blobs[source] = sb

        return source_blobs
    
    def get_constrained_catalogue(self) -> DataFrame:
        """
        Get the constrained catalogue containing only objects within the cutout area.
        
        :return: Constrained catalogue DataFrame
        :rtype: DataFrame
        """
        return self.constrained_catalogue

    def _constrain_catalogue(
            self,
            catalogue: DataFrame,
            size_arcmin: float,
            cutout_ra: float,
            cutout_dec: float,
            ra_col: str,
            dec_col: str,
        ):
        """
        Constrain the catalogue to only include objects within the cutout area.
        Computes a bounding box around the cutout center based on the size in arcminutes.

        :param catalogue: Pandas DataFrame containing the source catalogue
        :type catalogue: DataFrame
        :param size_arcmin: Size of the cutout in arcminutes
        :type size_arcmin: float
        :param cutout_ra: Right Ascension of the cutout center in degrees
        :type cutout_ra: float
        :param cutout_dec: Declination of the cutout center in degrees
        :type cutout_dec: float
        :param ra_col: Name of the RA column in the catalogue
        :type ra_col: str
        :param dec_col: Name of the Dec column in the catalogue
        :type dec_col: str
        :return: Constrained catalogue DataFrame
        :rtype: DataFrame
        """
        max_sep_arcsec = size_arcmin * 60 / 2

        cos_dec = np.cos(np.deg2rad(cutout_dec))
        if abs(cos_dec) < 1e-8:
            cos_dec = 1e-8

        delta_ra = max_sep_arcsec / cos_dec / 3600  # in degrees
        delta_dec = max_sep_arcsec / 3600  # in degrees

        ra_min = cutout_ra - delta_ra
        ra_max = cutout_ra + delta_ra
        dec_min = cutout_dec - delta_dec
        dec_max = cutout_dec + delta_dec

        ra_sorted = self._cache["ra_sorted"]
        dec_sorted = self._cache["dec_sorted"]
        sort_idx = self._cache["sort_idx"]

        left_idx = np.searchsorted(ra_sorted, ra_min, side="left")
        right_idx = np.searchsorted(ra_sorted, ra_max, side="right")

        if left_idx >= right_idx:
            return catalogue.iloc[0:0]

        candidate_sorted_idx = sort_idx[left_idx:right_idx]
        candidate_dec = dec_sorted[left_idx:right_idx]
        dec_mask = (candidate_dec >= dec_min) & (candidate_dec <= dec_max)

        if not np.any(dec_mask):
            return catalogue.iloc[0:0]

        filtered_idx = candidate_sorted_idx[dec_mask]
        return catalogue.iloc[filtered_idx]
    