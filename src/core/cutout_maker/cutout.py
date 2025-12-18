import os

from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.visualization import (AsinhStretch, ImageNormalize, PercentileInterval)
import dotenv
import numpy as np
import matplotlib.pyplot as plt
import gc

from ..mosaic import Mosaic


dotenv.load_dotenv()
RMS = float(os.getenv("RMS", 0.1))  # mJy/beam


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
        output_path = os.path.join(output_path, filename)
        
        # Save to FITS
        hdu = fits.PrimaryHDU(data=cutout.data, header=cutout_header)
        hdu.writeto(output_path, overwrite=True)

    def show_cutout(self, title: str, contour_levels: list = None, cmap='inferno', colorbar=False, save_path=None, save=False) -> None:
        """
        Show the cutout using matplotlib.

        :param title: Title for the plot
        :param contour_levels: List of contour levels in multiples of RMS noise (e.g., [3, 5, 10])
        :param cmap: Colormap for the image
        :param colorbar: Whether to show a colorbar
        :param save_path: Path to save the plot image (if save is True)
        :param save: Whether to save the plot instead of showing it
        """
        cutout = self.cutout
        data = cutout.data

        contour_levels = self._make_contour_levels(contour_levels) if contour_levels is not None else None

        norm = ImageNormalize(data, interval=PercentileInterval(99.5), stretch=AsinhStretch())
        fig = plt.figure(figsize=(10, 10))
        ax = plt.subplot(projection=cutout.wcs)
        im = ax.imshow(data, cmap=cmap, origin='lower', norm=norm)
        if colorbar:
            plt.colorbar(im, ax=ax, label='Intensity')
        if contour_levels is not None:
            ax.contour(data, levels=contour_levels, colors='white', linewidths=0.5)
        ax.set_title(title)
        ax.set_xlabel('RA (J2000)')
        ax.set_ylabel('Dec (J2000)')
        if save:
            if save_path is None:
                save_path = f"{title.replace(' ', '_')}_cutout.png"
            plt.savefig(save_path)
        else:
            plt.show()
        plt.close()

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
    
        self._cutout = Cutout2D(data, position, (self.size_pixels, self.size_pixels), wcs=wcs)
        return self._cutout

    def _make_contour_levels(self, levels: list) -> np.ndarray:
        """
        Generates contour levels based on RMS noise.
        
        :param levels: List of contour levels in multiples of RMS noise
        :return: List of contour levels in Jy/beam
        """
        return np.array(levels) * RMS * 1e-3  # Convert mJy/beam to Jy/beam
    
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
    