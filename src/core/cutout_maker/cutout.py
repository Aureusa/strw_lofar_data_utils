import os

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import dotenv
import numpy as np

import gc

from ..mosaic import Mosaic


class Cutout:
    """
    Class representing a cutout from a LOFAR mosaic.
    """
    def __init__(self, mosaic: Mosaic, ra: float, dec: float, size_arcmin: float = None, size_pixels: int= None):
        """
        Initialize a Cutout object.
        
        :param mosaic: Mosaic object from which to create the cutout
        :param ra: Right Ascension of the cutout center in degrees
        :param dec: Declination of the cutout center in degrees
        :param size_arcmin: Size of the cutout in arcminutes (optional) - should provide either this or size_pixels
        :param size_pixels: Size of the cutout in pixels (optional) - should provide either this or size_arcmin
        """
        if size_arcmin is None and size_pixels is None:
            raise ValueError("Either size_arcmin or size_pixels must be provided.")
        
        self.mosaic = mosaic
        self.ra = ra
        self.dec = dec

        self.size_arcmin = size_arcmin if size_arcmin is not None else self._calculate_size_arcmin(size_pixels)
        self.size_pixels = size_pixels if size_pixels is not None else self._calculate_size_pixels(size_arcmin)

        self._cutout = None

        # This flag indicates whether the cutout is valid (non-zero size)
        # This is temporally fix since sometimes Cutout2D returns zero-size
        # cutouts without raising an error - need to investigate further
        # The pipeline should check this flag after creating the cutout
        # and handle invalid cutouts appropriately. God forbit we have
        # to deal with this in multiple places.
        self.valid = self._check_valid_data_region()

    def _check_valid_data_region(self) -> bool:
        """
        Check if cutout center is within the valid (non-NaN) circular region.
        Uses angular separation without loading data.
        """        
        # Calculate angular separation between cutout center and mosaic center
        cutout_coord = SkyCoord(ra=self.ra*u.deg, dec=self.dec*u.deg, frame='icrs')
        mosaic_coord = SkyCoord(ra=self.mosaic.ra*u.deg, dec=self.mosaic.dec*u.deg, frame='icrs')
        
        separation = cutout_coord.separation(mosaic_coord).deg
        
        # Check if within valid radius (with some margin for cutout size)
        cutout_half_size_deg = self.size_arcmin / 60 / 2  # Convert arcmin to degrees, then half
        
        # Cutout is valid if its center + half its size is within the valid radius
        return (separation + cutout_half_size_deg) <= self.mosaic.valid_data_radius_deg

    def get_data(self) -> np.ndarray:
        """
        Get the cutout data as a 2D numpy array.
        """
        return self.cutout.data

    def get_header(self) -> fits.Header:
        """
        Get the cutout header as a FITS header.
        """
        return self.mosaic.load_header()
    
    def get_wcs(self) -> WCS:
        """
        Get the cutout WCS object.
        """
        return self.cutout.wcs

    @property
    def cutout(self) -> Cutout2D:
        """
        Get the Cutout2D object.
        """
        return self._create_cutout()

    def save_cutout(self, output_path: str, filename: str = None) -> None:
        """
        Save the cutout to a FITS file.

        :param output_path: Path to save the cutout FITS file
        :param filename: Optional filename for the cutout FITS file. Usually generated if not provided.
        """
        cutout = self.cutout

        # Create new header for cutout
        cutout_header = cutout.wcs.to_header()

        filename = self._generate_filename() if filename is None else filename

        # Make sure output path exists
        if not os.path.exists(output_path):
            os.makedirs(output_path)
            
        output_path = os.path.join(output_path, filename)
        
        # Save to FITS
        hdu = fits.PrimaryHDU(data=cutout.data, header=cutout_header)
        hdu.writeto(output_path, overwrite=True)

    def offload_data(self) -> None:
        """
        Offload the cutout data to save memory.
        """
        self._cutout = None
        self.mosaic.offload_data()
        gc.collect()

    def _create_cutout(self) -> Cutout2D:
        """
        Create and return the Cutout2D object.

        :return: Cutout2D object
        """
        if self._cutout is not None:
            return self._cutout
        
        # Load mosaic data and header
        data = self.mosaic.load_data()
        header = self.mosaic.load_header()
        wcs = WCS(header)

        # Offload the Mosaic data to save memory
        self.mosaic.offload_data()

        # Create cutout
        position = SkyCoord(ra=self.ra*u.deg, dec=self.dec*u.deg, frame='icrs')
    
        self._cutout = Cutout2D(
            data,
            position,
            (self.size_pixels, self.size_pixels),
            wcs=wcs
        )
        return self._cutout
    
    def _generate_filename(self) -> str:
        """
        Generate a filename for the cutout FITS file.
        """
        ra_str = f"{self.ra:.5f}".replace('.', 'p')
        dec_str = f"{self.dec:.5f}".replace('.', 'p').replace('-', 'm')
        return f"{self.mosaic.field_name}_cutout_RA{ra_str}_Dec{dec_str}_Size{self.size_pixels}px.fits"  

    def _calculate_size_arcmin(self, size_pixels: int) -> float:
        """
        Calculate size in arcminutes from size in pixels.
        
        :param size_pixels: Size in pixels
        :return: Size in arcminutes
        """
        pixel_scale_deg = abs(self.mosaic.header['CDELT1'])  # degrees per pixel
        size_deg = size_pixels * pixel_scale_deg
        return size_deg * 60  # convert to arcminutes
    
    def _calculate_size_pixels(self, size_arcmin: float) -> int:
        """
        Calculate size in pixels from size in arcminutes.
        
        :param size_arcmin: Size in arcminutes
        :return: Size in pixels
        """
        pixel_scale_deg = abs(self.mosaic.header['CDELT1'])  # degrees per pixel
        size_deg = size_arcmin / 60  # convert to degrees
        return int(size_deg / pixel_scale_deg)
    