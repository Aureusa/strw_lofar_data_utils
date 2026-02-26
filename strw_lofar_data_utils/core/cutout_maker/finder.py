from ..mosaic import Mosaic


def find_candidate_mosaics(ra: float, dec: float, list_of_mosaics: list[Mosaic]) -> list[Mosaic]:
    """
    Find candidate mosaics that cover the given RA and Dec.

    :param ra: Right Ascension in degrees
    :param dec: Declination in degrees
    :param list_of_mosaics: List of Mosaic objects
    :return: List of Mosaic objects that cover the position, or empty list if not found
    """
    candidate_mosaics = []
    for mosaic in list_of_mosaics:
        if mosaic.is_in_coverage(ra, dec):
            candidate_mosaics.append(mosaic)
    return candidate_mosaics

def chose_from_candidate_mosaics(ra: float, dec: float, candidate_mosaics: list[Mosaic]) -> Mosaic:
    """
    Choose the best mosaic from a list of candidate mosaics covering the given RA and Dec.
    The best mosaic is chosen based on the smallest separation from the center.

    :param ra: Right Ascension in degrees
    :param dec: Declination in degrees
    :param candidate_mosaics: List of candidate Mosaic objects
    :return: The chosen Mosaic object
    """
    separations = []
    for mosaic in candidate_mosaics:
        separations.append(
            mosaic.find_seperation_from_center(ra, dec)
        )

    # Chose the mosaic with the smallest separation from the center
    best_mosaic = None
    if separations:
        best_index = separations.index(min(separations))
        best_mosaic = candidate_mosaics[best_index]

    return best_mosaic
