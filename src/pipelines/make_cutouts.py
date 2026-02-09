from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm

from ..core.cutout_maker import make_cutout, find_candidate_mosaics, make_value_added_cutout, chose_from_candidate_mosaics
from ..core.mosaic import get_list_of_mosaics


def _process_single_cutout(ra_dec, mosaics, size_arcmin, size_pixels, data_folder, save, value_added):
    """Helper function to process a single cutout (for parallel execution)."""
    ra, dec = ra_dec
    candidate_mosaic = find_candidate_mosaics(ra, dec, mosaics)
    mosaic = chose_from_candidate_mosaics(ra, dec, candidate_mosaic)

    
    if mosaic is None:
        return None, f"No mosaic found covering RA: {ra}, Dec: {dec}"
    
    if value_added:
        cutout = make_value_added_cutout(
            mosaic,
            ra,
            dec,
            size_arcmin=size_arcmin,
            size_pixels=size_pixels,
        )
    else:
        cutout = make_cutout(
            mosaic,
            ra,
            dec,
            size_arcmin=size_arcmin,
            size_pixels=size_pixels,
        )

    if cutout is None:
        return None, f"Cutout at RA: {ra}, Dec: {dec} is outside valid data region."
    
    if save:
        cutout.save_cutout(output_path=data_folder)
    
    return cutout, None


def generate_cutouts(
        ra_dec_list: list[tuple[float, float]],
        size_arcmin: float = None,
        size_pixels: int = None,
        data_folder: str = None,
        save: bool = False,
        mosaic_coverage_file: str = 'default',
        value_added: bool = False,
        n_workers: int = None,
        verbose: bool = False
    ) -> list:
    """
    Generate cutouts for a list of RA and Dec positions.

    :param ra_dec_list: List of tuples containing (RA, Dec) in degrees
    :param size_arcmin: Size of the cutout in arcminutes (optional)
    :param size_pixels: Size of the cutout in pixels (optional)
    :param data_folder: Folder to save cutouts if save is True
    :param save: Whether to save the cutouts to disk
    :param mosaic_coverage_file: Path to the mosaic coverage CSV file. If 'default'
    uses the coverage file in the `/path/to/repo/data/mosaic_coverage/lotss_dr2_mosaic_coverage.csv`.
    In general I don't see a reason to change this.
    :param value_added: Whether to create value-added cutouts
    :param n_workers: Number of parallel workers (None = use all CPUs)
    :return: List of Cutout objects
    """
    cutouts = []
    mosaics = get_list_of_mosaics(mosaic_coverage_file)
    
    # Create partial function with fixed parameters
    process_func = partial(
        _process_single_cutout,
        mosaics=mosaics,
        size_arcmin=size_arcmin,
        size_pixels=size_pixels,
        data_folder=data_folder,
        save=save,
        value_added=value_added,
    )
    
    # Parallel processing
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_func, ra_dec): ra_dec for ra_dec in ra_dec_list}
        
        for future in tqdm(as_completed(futures), total=len(ra_dec_list), desc="Generating cutouts"):
            cutout, error = future.result()
            if error:
                if verbose:
                    print(error)
            elif cutout is not None:
                cutouts.append(cutout)
    
    return cutouts
