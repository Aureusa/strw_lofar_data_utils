import os
import gc

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

from .utils import find_field_path


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
            dec_size: float,
            release: str = "DR2"
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
        :param ra_size: Size of the mosaic in RA direction in degrees
        :param dec_size: Size of the mosaic in Dec direction in degrees
        :param release: Data release ("DR2" or "DR3")
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
        self.release = release

        self._header = None
        self._data = None

        # Find the mosaic paths
        self._field_path = find_field_path(self.field_name, self.release)
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
        blanked_path = os.path.join(
            self._field_path,
            self.field_name,
            mosaic_blanked
        )
        if not os.path.exists(blanked_path):
            raise FileNotFoundError(f"Blanked mosaic FITS file not found at expected path: {blanked_path}")
        return blanked_path
    
    def _get_center_coords(self) -> SkyCoord:
        """
        Get the SkyCoord of the mosaic center.

        :return: SkyCoord of the mosaic center
        """
        return SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg)
