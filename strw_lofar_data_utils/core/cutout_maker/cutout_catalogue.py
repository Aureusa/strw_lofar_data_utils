from pandas.core.frame import DataFrame
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

from .source_blob import SourceBlob
from .cutout import Cutout


class CutoutCatalogue:
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

        unique_sources = self.constrained_catalogue[target_column].unique()

        # Create SourceBlob instances for each unique source
        source_blobs = {}
        for source in unique_sources:
            source_data = self.constrained_catalogue[self.constrained_catalogue[target_column] == source]
            ra_list = source_data[self.ra_col].tolist()
            dec_list = source_data[self.dec_col].tolist()
            sb = SourceBlob(ra_list, dec_list)
            in_cutout = sb.check_if_in_cutout(self.cutout_wcs, self.cutout_shape)
            
            # If all are in cutout, add to dict
            if all(in_cutout):
                source_blobs[source] = sb
            elif any(in_cutout):
                # If some are in cutout, create new SourceBlob with only those
                ra_in = [ra for ra, inside in zip(ra_list, in_cutout) if inside]
                dec_in = [dec for dec, inside in zip(dec_list, in_cutout) if inside]
                sb_in = SourceBlob(ra_in, dec_in)
                sb_in.check_if_in_cutout(self.cutout_wcs, self.cutout_shape)  # Should all be True now
                source_blobs[source] = sb_in
            else:
                # None are in cutout, skip
                continue

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

        # Create a box around the given RA, Dec
        center_coord = SkyCoord(ra=cutout_ra * u.deg, dec=cutout_dec * u.deg, frame='icrs')
        delta_ra = max_sep_arcsec / np.cos(np.deg2rad(cutout_dec)) / 3600  # in degrees
        delta_dec = max_sep_arcsec / 3600  # in degrees

        # Extract maximum and minimum RA and Dec
        ra_min = center_coord.ra.deg - delta_ra
        ra_max = center_coord.ra.deg + delta_ra
        dec_min = center_coord.dec.deg - delta_dec
        dec_max = center_coord.dec.deg + delta_dec

        constrained_catalogue = catalogue[
            (catalogue[ra_col] >= ra_min) & (catalogue[ra_col] <= ra_max) &
            (catalogue[dec_col] >= dec_min) & (catalogue[dec_col] <= dec_max)
        ]
        return constrained_catalogue
    