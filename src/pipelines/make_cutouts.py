from tqdm import tqdm

from ..core.cutout_maker import make_cutout, find_mosaic
from ..core.mosaic import get_list_of_mosaics

def generate_cutouts(
        ra_dec_list: list[tuple[float, float]],
        size_arcmin: float = None,
        size_pixels: int = None,
        data_folder: str = None,
        save: bool = False
    ) -> list:
    """
    Generate cutouts for a list of RA and Dec positions.

    :param ra_dec_list: List of tuples containing (RA, Dec) in degrees
    :param size_arcmin: Size of the cutout in arcminutes (optional)
    :param size_pixels: Size of the cutout in pixels (optional)
    :param data_folder: Folder to save cutouts if save is True
    :param save: Whether to save the cutouts to disk
    :return: List of Cutout objects
    """
    cutouts = []
    mosaics = get_list_of_mosaics()

    for ra, dec in tqdm(ra_dec_list, desc="Generating cutouts"):
        mosaic = find_mosaic(ra, dec, mosaics)
        if mosaic is not None:
            cutout = make_cutout(mosaic, ra, dec, size_arcmin=size_arcmin, size_pixels=size_pixels)

            if save:
                cutout.save_cutout(output_path=data_folder)
                cutout.offload_data()

            cutouts.append(cutout)
        else:
            print(f"No mosaic found covering RA: {ra}, Dec: {dec}")

    return cutouts
