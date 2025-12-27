from ..mosaic import Mosaic
from .cutout import Cutout
from .value_added_cutout import ValueAddedCutout
from .astro_object import AstroObject

def make_cutout(mosaic: Mosaic, ra: float, dec: float, size_arcmin: float = None, size_pixels: int = None):
    """
    Create a Cutout object for the given mosaic at specified RA, Dec.

    :param mosaic: Mosaic object
    :param ra: Right Ascension in degrees
    :param dec: Declination in degrees
    :param size_arcmin: Size of the cutout in arcminutes (optional)
    :param size_pixels: Size of the cutout in pixels (optional)
    :return: Cutout object
    """
    cutout = Cutout(mosaic, ra, dec, size_arcmin=size_arcmin, size_pixels=size_pixels)
    return cutout

def make_value_added_cutout(
    mosaic: Mosaic,
    ra: float,
    dec: float,
    size_arcmin: float = None,
    size_pixels: int = None,
    ):
    """
    Create a ValueAddedCutout object for the given mosaic at specified RA, Dec.

    :param mosaic: Mosaic object
    :param ra: Right Ascension in degrees
    :param dec: Declination in degrees
    :param size_arcmin: Size of the cutout in arcminutes (optional)
    :param size_pixels: Size of the cutout in pixels (optional)
    :return: ValueAddedCutout object
    """
    cutout = ValueAddedCutout(
        mosaic,
        ra,
        dec,
        size_arcmin=size_arcmin,
        size_pixels=size_pixels,
    )
    return cutout