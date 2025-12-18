import os
from astropy.io import fits
import argparse
import pandas as pd
import yaml
from tqdm import tqdm
import dotenv

dotenv.load_dotenv()

# Get environment variables
BASE_DIR = os.getenv("BASE_DIR", "/disks/paradata/shimwell/LoTSS-DR2/mosaics")
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
    

def print_config_summary(
        base_dir: str, field_list: list, save_csv: bool, save_path: str, verbose: bool
    ) -> None:
    """
    Print a summary of the configuration settings.

    :param config: Configuration dictionary
    :param save_csv: Whether to save the output DataFrame to a CSV file
    :param save_path: Path to save the CSV file
    :param verbose: Whether to print verbose output
    """
    print("Configuration Summary:")
    print(f"Base Directory: {base_dir}")
    print(f"Fields: {field_list}")
    print(f"Save CSV: {save_csv}")
    print(f"CSV Save Path: {save_path}")
    print(f"Verbose: {verbose}")
    

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
        verbose: bool = True
    ) -> None:
    """
    Main function to extract mosaic coverage information and save to CSV.
    1. Iterate over specified fields and pointings.
    2. Extract coverage info from each mosaic FITS file.
    3. Compile results into a DataFrame and save as CSV if required.
    4. Print verbose output if enabled.
    
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
    :param verbose: Whether to print verbose output. Default is True.
    """
    # Collect all mosaic info
    all_mosaics = []
    for ra_field in field_list:
        field_dir = os.path.join(base_dir, ra_field)
        if os.path.exists(field_dir):
            for pointing in tqdm(os.listdir(field_dir), desc=f"Processing {ra_field}") if verbose else os.listdir(field_dir):
                mosaic_path = os.path.join(field_dir, pointing, "mosaic-blanked.fits")
                if os.path.exists(mosaic_path):
                    info = get_mosaic_coverage(mosaic_path)
                    all_mosaics.append(info)

    df = pd.DataFrame(all_mosaics)

    if verbose:
        print(df.head())

    if save_csv:
        df.to_csv(save_path, index=False)
        print(f"\nSaved {len(df)} mosaics to {save_path}")

if __name__ == "__main__":
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

    print_config_summary(BASE_DIR, FIELD_LIST, save_csv, save_path, verbose)

    main(
        BASE_DIR,
        FIELD_LIST,
        save_csv=save_csv,
        save_path=save_path,
        verbose=verbose
    )
    