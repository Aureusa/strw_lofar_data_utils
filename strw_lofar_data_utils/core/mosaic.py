import os
import glob
import dotenv
import gc

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd
from tqdm import tqdm
from .utils import pixel_2_arcmin_size, arcmin_2_pixel_size


dotenv.load_dotenv()
BASE_DIR = os.getenv("BASE_DIR", "/disks/paradata/shimwell/LoTSS-DR2/mosaics")
RA0h_field = os.getenv("RA0h_field", "RA0h_field")
RA13h_field = os.getenv("RA13h_field", "RA13h_field")
FIELD_LIST = [RA0h_field, RA13h_field]


class Mosaic:
    """
    Class representing a LOFAR mosaic.
    """
    def __init__(
            self,
            field_name: str,
            ra: float, dec: float,
            ra_min: float, 
            ra_max: float,
            dec_min: float,
            dec_max: float,
            ra_size: float,
            dec_size: float
        ):
        """
        Initialize a Mosaic object.
        
        :param field_name: Name of the mosaic field
        :param ra: Central RA of the mosaic in degrees
        :param dec: Central Dec of the mosaic in degrees
        :param ra_min: Minimum RA of the mosaic coverage in degrees
        :param ra_max: Maximum RA of the mosaic coverage in degrees
        :param dec_min: Minimum Dec of the mosaic coverage in degrees
        :param dec_max: Maximum Dec of the mosaic coverage in degrees
        :param ra_size: Size of the mosaic in RA direction in degrself.raees
        :param dec_size: Size of the mosaic in Dec direction in degrees
        """
        self.field_name = field_name
        self.ra = ra
        self.dec = dec
        self.ra_min = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max
        self.ra_size = ra_size
        self.dec_size = dec_size

        self._header = None
        self._data = None

        # Find the mosaic paths
        self._field_path = self._find_field_path()
        self._mosaic_blanked_path = self._find_candidate_mosaics_blanked_path()

        # Define the valid data region (non-NAN) which is a circle around the center
        # with radius equal to half the smaller of ra_size and dec_size
        self.valid_data_radius_deg = min(self.ra_size, self.dec_size) / 2

    @property
    def wcs(self):
        """
        Get the WCS object for the mosaic.
        """
        if self._header is None:
            self.load_header()
        return WCS(self._header)

    @property
    def header(self):
        """
        Get the FITS header for the mosaic.
        """
        if self._header is None:
            self.load_header()
        return self._header
    
    @property
    def data(self):
        """
        Get the FITS data for the mosaic.
        Must call load_data() first.
        """
        if self._data is None:
            print("Data not loaded. Call load_data() to load it first.")
        return self._data
    
    def find_empty_regions(
            self,
            pybdsf_catalog: pd.DataFrame,
            ra_key: str,
            dec_key: str,
            size_arcmin: float = None,
            size_pixels: int = None,
            max_sources: int = 3,
            grid_step_arcmin: float = None,
            grid_step_pixels: int = None,
            min_separation_arcmin: float = None
        ) -> pd.DataFrame:
        """
        Find empty regions in the mosaic with few or no sources from the catalog.
        
        :param pybdsf_catalog: DataFrame containing the PyBDSF catalog
        :param ra_key: Column name for RA in the catalog
        :param dec_key: Column name for Dec in the catalog
        :param size_arcmin: Size of the region to check in arcminutes
        :param size_pixels: Size of the region to check in pixels
        :param max_sources: Maximum number of sources allowed in an "empty" region (default: 0)
        :param grid_step_arcmin: Step size for grid search in arcminutes (default: size_arcmin/2)
        :param grid_step_pixels: Step size for grid search in pixels (default: size_pixels/2)
        :param min_separation_arcmin: Minimum separation between empty regions (default: size_arcmin)
        :return: DataFrame with columns: ra, dec, num_sources
        """
        if size_arcmin is None and size_pixels is None:
            raise ValueError("Either size_arcmin or size_pixels must be provided.")
        if size_pixels is None:
            size_pixels = arcmin_2_pixel_size(size_arcmin, self.header)
        if size_arcmin is None:
            size_arcmin = pixel_2_arcmin_size(size_pixels, self.header)
        
        # Set default grid step to half the region size
        if grid_step_arcmin is None and grid_step_pixels is None:
            grid_step_arcmin = size_arcmin / 2
        if grid_step_pixels is None:
            grid_step_pixels = arcmin_2_pixel_size(grid_step_arcmin, self.header)
        if grid_step_arcmin is None:
            grid_step_arcmin = pixel_2_arcmin_size(grid_step_pixels, self.header)
        
        # Set default minimum separation
        if min_separation_arcmin is None:
            min_separation_arcmin = size_arcmin

        # Convert catalog to SkyCoord
        catalog_coords = SkyCoord(
            ra=pybdsf_catalog[ra_key].values * u.deg,
            dec=pybdsf_catalog[dec_key].values * u.deg
        )

        # Create a WCS object for coordinate transformations
        wcs = WCS(self.header)
        
        # Calculate search radius (half diagonal of the box)
        search_radius_arcmin = size_arcmin * np.sqrt(2) / 2
        search_radius = search_radius_arcmin * u.arcmin

        # Create grid of test positions
        ra_grid = np.arange(self.ra_min, self.ra_max, grid_step_arcmin / 60)
        dec_grid = np.arange(self.dec_min, self.dec_max, grid_step_arcmin / 60)
        
        empty_regions = []
        
        for ra in tqdm(ra_grid, desc="Scanning RA-DEC grid"):
            for dec in dec_grid:
                # Check if the center is within valid data region
                center = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
                separation_from_mosaic_center = center.separation(
                    SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg)
                )
                
                if separation_from_mosaic_center.deg > self.valid_data_radius_deg:
                    continue
                
                # Find sources within the search radius
                separations = center.separation(catalog_coords)
                sources_in_region = np.sum(separations < search_radius)
                
                # If region is empty enough, add it to the list
                if sources_in_region <= max_sources:
                    empty_regions.append({
                        'ra': ra,
                        'dec': dec,
                        'num_sources': sources_in_region
                    })
        
        # Convert to DataFrame
        empty_regions_df = pd.DataFrame(empty_regions)
        
        if len(empty_regions_df) == 0:
            return empty_regions_df
        
        # Filter out regions that are too close to each other
        # Keep the ones with fewer sources
        empty_regions_df = empty_regions_df.sort_values('num_sources')
        filtered_regions = []
        
        for _, row in tqdm(empty_regions_df.iterrows(), total=len(empty_regions_df), desc="Filtering close regions"):
            current_coord = SkyCoord(ra=row['ra'] * u.deg, dec=row['dec'] * u.deg)
            
            # Check if this region is far enough from already selected regions
            too_close = False
            for selected in filtered_regions:
                selected_coord = SkyCoord(ra=selected['ra'] * u.deg, dec=selected['dec'] * u.deg)
                if current_coord.separation(selected_coord).arcmin < min_separation_arcmin:
                    too_close = True
                    break
            
            if not too_close:
                filtered_regions.append(row.to_dict())
        
        return pd.DataFrame(filtered_regions)
    
    def find_seperation_from_center(self, ra: float, dec: float) -> float:
        """
        Find the angular separation from the mosaic center to the given RA and Dec.

        :param ra: Right Ascension in degrees
        :param dec: Declination in degrees
        :return: Separation in degrees
        """
        center_coord = SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg)
        target_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        separation = center_coord.separation(target_coord)
        return separation.deg

    def is_in_coverage(self, ra: float, dec: float) -> bool:
        """
        Check if the given RA and Dec are within the mosaic coverage.

        :param ra: Right Ascension in degrees
        :param dec: Declination in degrees
        :return: True if within coverage, False otherwise
        """
        return (self.ra_min <= ra <= self.ra_max) and (self.dec_min <= dec <= self.dec_max)

    def convert_pixels_to_world(self, x_pixel: list[float], y_pixel: list[float]) -> list[tuple[float, float]]:
        """
        Convert pixel coordinates to world coordinates (RA, Dec).

        :param x_pixel: X pixel coordinate
        :type x_pixel: list of floats
        :param y_pixel: Y pixel coordinate
        :type y_pixel: list of floats
        :return: List of tuples of (RA, Dec) in degrees
        """
        wcs = WCS(self.header)

        x_array = np.asarray(x_pixel, dtype=np.float64)
        y_array = np.asarray(y_pixel, dtype=np.float64)
        ra_array, dec_array = wcs.wcs_pix2world(x_array, y_array, 0)
        return list(zip(ra_array.tolist(), dec_array.tolist()))
    
    def load_header(self):
        """
        Load the FITS header for the mosaic.

        :return: FITS header
        """
        if self._header is None:
            fits_path = self._mosaic_blanked_path
            with fits.open(fits_path) as hdul:
                self._header = hdul[0].header
        return self._header
    
    def load_data(self):
        """
        Load the FITS data for the mosaic.

        :return: FITS data
        """
        if self._data is None:
            fits_path = self._mosaic_blanked_path
            with fits.open(fits_path) as hdul:
                self._data = hdul[0].data
        return self._data

    def get_data_shape(self):
        """
        Get the shape of the FITS data for the mosaic without loading the entire data into memory.
    
        :return: Shape of the FITS data
        """
        if self._data is not None:
            return self._data.shape
        
        fits_path = self._mosaic_blanked_path
        with fits.open(fits_path) as hdul:
            data_shape = hdul[0].data.shape
        return data_shape
    
    def offload_data(self):
        """
        Offload the FITS data from memory.
        """
        self._data = None
        self._header = None
        gc.collect()

    def _find_candidate_mosaics_blanked_path(self):
        """
        Find the path to the blanked mosaic FITS file.
        Assumes a standard naming convention - 'mosaic-blanked.fits'.

        :return: Path to the blanked mosaic FITS file
        """
        mosaic_blanked = "mosaic-blanked.fits"
        return os.path.join(
            self._field_path,
            self.field_name,
            mosaic_blanked
        )
    
    def _get_center_coords(self) -> SkyCoord:
        """
        Get the SkyCoord of the mosaic center.

        :return: SkyCoord of the mosaic center
        """
        return SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg)

    def _find_field_path(self) -> str:
        """
        Find the base path for the mosaic field.

        :return: Base path for the mosaic field
        """
        # Get all subdirectories in BASE_DIR + the two field names
        search_paths = [
            os.path.join(BASE_DIR, RA0h_field),
            os.path.join(BASE_DIR, RA13h_field)
        ]

        for base_path in search_paths:
            subdirs = glob.glob(os.path.join(base_path, "*"))
            for subdir in subdirs:
                if os.path.basename(subdir) == self.field_name:
                    return base_path
        raise FileNotFoundError(f"Field {self.field_name} not found in expected directories.")
    

def get_list_of_mosaics(mosaic_coverage_file: str = 'default') -> list[Mosaic]:
    """
    Get the list of Mosaic objects from the mosaic coverage CSV file. It assumes
    that the CSV file has columns: field_name, ra, dec, ra_min, ra_max, dec_min, dec_max, ra_size, dec_size.
    Usually this file is generated by running the `get_mosaic_coverage.sh` script in /scripts folder.
    
    :param mosaic_coverage_file: Path to the mosaic coverage CSV file
    :return: List of Mosaic objects
    """
    if mosaic_coverage_file == 'default':
        mosaic_coverage_file = os.path.join(
            os.path.dirname(__file__),
            '..',
            '..',
            'data',
            'mosaic_coverage',
            'lotss_dr2_mosaic_coverage.csv'
        )
    df = pd.read_csv(mosaic_coverage_file)

    mosaics = []
    for _, row in df.iterrows():
        mosaic = Mosaic(
            field_name=row['field_name'],
            ra=row['ra'],
            dec=row['dec'],
            ra_min=row['ra_min'],
            ra_max=row['ra_max'],
            dec_min=row['dec_min'],
            dec_max=row['dec_max'],
            ra_size=row['ra_size'],
            dec_size=row['dec_size']
        )
        mosaics.append(mosaic)
    return mosaics


def get_mosaic_by_field_name(field_name: str, mosaic_coverage_file: str = 'default') -> Mosaic:
    """
    Get a Mosaic object by its field name.

    :param field_name: Name of the mosaic field
    :param mosaic_coverage_file: Path to the mosaic coverage CSV file
    :return: Mosaic object with the given field name
    """
    if mosaic_coverage_file == 'default':
        mosaic_coverage_file = os.path.join(
            os.path.dirname(__file__),
            '..',
            '..',
            'data',
            'mosaic_coverage',
            'lotss_dr2_mosaic_coverage.csv'
        )
    df = pd.read_csv(mosaic_coverage_file)

    # Find the row with the given field name
    row = df[df['field_name'] == field_name]
    if row.empty:
        raise ValueError(f"Mosaic with field name {field_name} not found in coverage file.")
    row = row.iloc[0]
    mosaic = Mosaic(
        field_name=row['field_name'],
        ra=row['ra'],
        dec=row['dec'],
        ra_min=row['ra_min'],
        ra_max=row['ra_max'],
        dec_min=row['dec_min'],
        dec_max=row['dec_max'],
        ra_size=row['ra_size'],
        dec_size=row['dec_size']
    )
    return mosaic
    