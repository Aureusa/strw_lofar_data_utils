import os
from astropy.io import fits
import argparse
import pandas as pd
from pathlib import Path
import yaml
from tqdm import tqdm
import dotenv
from logging import Logger

from ...logger import setup_logging

dotenv.load_dotenv(dotenv_path=Path(__file__).resolve().parents[3] / ".env")

# Get environment variables
DR2_BASE_DIR = os.getenv("DR2_BASE_DIR", "/disks/paradata/shimwell/LoTSS-DR2/mosaics")
DR3_BASE_DIR = os.getenv("DR3_BASE_DIR", "/disks/paradata/shimwell/Beyond-DR2/mosaics/LoTSS-DR3-mosaics/")
RA0h_field = os.getenv("RA0h_field", "RA0h_field")
RA13h_field = os.getenv("RA13h_field", "RA13h_field")
FIELD_LIST = [RA0h_field, RA13h_field]


def load_config(config_path: str) -> dict:
    """
    Load configuration from a YAML file.

    :param config_path: Path to the YAML configuration file
    :return: Configuration dictionary
    """
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)
    

def log_config_summary(
        base_dir: str, field_list: list, save_csv: bool, save_path: str, logger: Logger = None, verbose: bool = True
    ) -> None:
    """
    Log a summary of the configuration settings.

    :param config: Configuration dictionary
    :param save_csv: Whether to save the output DataFrame to a CSV file
    :param save_path: Path to save the CSV file
    :param verbose: Whether to Log verbose output
    """
    if logger is None:
        logger = setup_logging(name="strw_lofar_data_utils.core.mosaic.get_mosaic_coverage")
    info = "Configuration Summary:"
    info += f"\nBase Directory: {base_dir}"
    info += f"\nFields: {field_list}"
    info += f"\nSave CSV: {save_csv}"
    info += f"\nCSV Save Path: {save_path}"
    info += f"\nVerbose: {verbose}"
    logger.info(info)

def get_mosaic_coverage(fits_path: str) -> dict:
    """
    Extract RA, Dec, and field size from FITS header
    
    :param fits_path: Path to the mosaic FITS file
    :return: Dictionary with field name, RA, Dec, RA/Dec min/max, and RA/Dec size
    """
    with fits.open(fits_path) as hdul:
        header = hdul[0].header
        
    ra_center = header['CRVAL1']  # RA in degrees
    dec_center = header['CRVAL2']  # Dec in degrees
    
    # Calculate field size
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']
    cdelt1 = abs(header['CDELT1'])  # pixel scale in degrees
    cdelt2 = abs(header['CDELT2'])
    
    ra_size = naxis1 * cdelt1  # field width in degrees
    dec_size = naxis2 * cdelt2  # field height in degrees
    
    # Calculate boundaries
    ra_min = ra_center - ra_size / 2
    ra_max = ra_center + ra_size / 2
    dec_min = dec_center - dec_size / 2
    dec_max = dec_center + dec_size / 2
    
    return {
        'field_name': header.get('OBJECT', 'Unknown').strip(),
        'ra': ra_center,
        'dec': dec_center,
        'ra_min': ra_min,
        'ra_max': ra_max,
        'dec_min': dec_min,
        'dec_max': dec_max,
        'ra_size': ra_size,
        'dec_size': dec_size
    }

def main(
        base_dir: str,
        field_list: list,
        save_csv: bool = True,
        save_path: str = 'lotss_dr2_mosaic_coverage.csv',
        logger=None,
        verbose: bool = True
    ) -> None:
    """
    Main function to extract mosaic coverage information and save to CSV.
    1. Iterate over specified fields and pointings.
    2. Extract coverage info from each mosaic FITS file.
    3. Compile results into a DataFrame and save as CSV if required.
    4. Log verbose output if enabled.
    
    :param base_dir: Base directory containing LOFAR DR2 mosaics. This
    should be set in the configuration YAML file. The default value is
    "/disks/paradata/shimwell/LoTSS-DR2/mosaics", which is valid for the
    STRW cluster.
    :param field_list: List of RA fields to process. This should be set
    in the configuration YAML file. Default fields are "RA0h_field" and
    "RA13h_field".
    :param save_csv: Whether to save the output DataFrame to a CSV file.
    Default is True.
    :param save_path: Path to save the CSV file. Default is
    'lotss_dr2_mosaic_coverage.csv'.
    :param verbose: Whether to log verbose output. Default is True.
    """
    # Collect all mosaic info
    all_mosaics = []
    for pointing in tqdm(os.listdir(base_dir), desc="Processing mosaics") if verbose else os.listdir(base_dir):
        mosaic_path = os.path.join(base_dir, pointing, "mosaic-blanked.fits")
        if os.path.exists(mosaic_path):
            info = get_mosaic_coverage(mosaic_path)
            all_mosaics.append(info)

    df = pd.DataFrame(all_mosaics)

    if verbose:
        logger.info(df.head())

    if save_csv:
        df.to_csv(save_path, index=False)
        logger.info(f"Saved {len(df)} mosaics to {save_path}")

if __name__ == "__main__":
    logger = setup_logging(name="strw_lofar_data_utils.core.mosaic.get_mosaic_coverage")

    parser = argparse.ArgumentParser(
        description="LOFAR Mosaic Coverage Extraction"
    )
    parser.add_argument(
        "--no-save-csv",
        action='store_false',
        dest='save_csv',
        help="Disable saving the output to CSV"
    )
    parser.add_argument(
        "--save-path",
        type=str,
        default='lotss_dr2_mosaic_coverage.csv',
        help="Path to save the CSV file (default: lotss_dr2_mosaic_coverage.csv)"
    )
    parser.add_argument(
        "--no-verbose",
        action='store_false',
        dest='verbose',
        help="Disable verbose output"
    )

    # Set default values
    parser._defaults['verbose'] = True
    parser._defaults['save_csv'] = True

    args = parser.parse_args()

    # Load configuration
    save_csv = args.save_csv
    save_path = args.save_path
    verbose = args.verbose

    log_config_summary(DR3_BASE_DIR, FIELD_LIST, save_csv, save_path, logger, verbose)

    main(
        DR3_BASE_DIR,
        FIELD_LIST,
        save_csv=save_csv,
        save_path=save_path,
        logger=logger,
        verbose=verbose
    )
    