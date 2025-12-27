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

from .cutout import Cutout
from .astro_object import AstroObject
from ..mosaic import Mosaic
from .astro_object import find_nearby_objects


dotenv.load_dotenv()


class ValueAddedCutout(Cutout):
    """
    Class representing a value-added cutout from a LOFAR mosaic.
    Extends the Cutout class with additional functionalities.
    """
    def __init__(self, mosaic: Mosaic, ra: float, dec: float, size_arcmin: float = None, size_pixels: int= None):
        """
        Initialize a ValueAddedCutout object.
        
        :param mosaic: Mosaic object from which to create the cutout
        :param ra: Right Ascension of the cutout center in degrees
        :param dec: Declination of the cutout center in degrees
        :param size_arcmin: Size of the cutout in arcminutes (optional) - should provide either this or size_pixels
        :param size_pixels: Size of the cutout in pixels (optional) - should provide either this or size_arcmin
        """
        super().__init__(mosaic, ra, dec, size_arcmin, size_pixels)
        header = self.mosaic.header
        self.value_added_data = find_nearby_objects(
            ra_dec=(ra, dec),
            header=header,
            catalogue_path=os.getenv('VALUE_ADDED_CATALOGUE_PATH'),
            size_arcmin=size_arcmin,
            size_pixels=size_pixels
        )

    def get_objects_in_cutout(self) -> dict[str, list[AstroObject]]:
        """
        Get the value-added objects that fall within the cutout area.

        :return: Dictionary with keys as object types and values as lists of AstroObject instances within the cutout
        """
        objects_in_cutout = {}
        cutout_wcs = self.cutout.wcs
        cutout_shape = self.cutout.data.shape

        for obj_name, astro_object in self.value_added_data.items():
            objects_in_cutout[obj_name] = []
            in_cutout = astro_object.check_if_in_cutout(cutout_wcs, cutout_shape)

            if not all(in_cutout):
                print(f"Warning: Some {obj_name} components are outside the cutout area.")
                objects_in_cutout[obj_name] = astro_object.get_pixel_positions()
            elif not any(in_cutout):
                objects_in_cutout[obj_name] = []
            else:
                objects_in_cutout[obj_name] = astro_object.get_pixel_positions()

        return objects_in_cutout

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

        objects_in_cutout = self.get_objects_in_cutout()

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

        # Plot objects in cutout as circles on the x y positions
        for obj_name, positions in objects_in_cutout.items():
            if positions:
                x_positions, y_positions = zip(*positions)
                ax.scatter(
                    x_positions,
                    y_positions,
                    s=100,
                    edgecolors='cyan',
                    facecolors='none',
                    marker='o',
                    label=obj_name
                )
        if objects_in_cutout:
            ax.legend(loc='upper right')
            
        if save:
            if save_path is None:
                save_path = f"{title.replace(' ', '_')}_cutout.png"
            plt.savefig(save_path)
        else:
            plt.show()
        plt.close()
    