from ..mosaic import Mosaic


def find_mosaic(ra: float, dec: float, list_of_mosaics: list[Mosaic]) -> Mosaic:
    """
    Find the mosaic that covers the given RA and Dec.

    :param ra: Right Ascension in degrees
    :param dec: Declination in degrees
    :param list_of_mosaics: List of Mosaic objects
    :return: Mosaic object that covers the position, or None if not found
    """
    for mosaic in list_of_mosaics:
        if mosaic.is_in_coverage(ra, dec):
            return mosaic
    return None
