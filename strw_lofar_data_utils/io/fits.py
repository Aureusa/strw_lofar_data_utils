import astropy.io.fits as fits


def read_fits_header(fits_path: str, idx: int = 0) -> fits.Header:
    """
    Read the header of a FITS file.
    
    :param fits_path: Path to the FITS file
    :param idx: Index of the HDU to read the header from (default is 0 for the primary HDU)
    :return: FITS header object
    """
    with fits.open(fits_path) as hdul:
        header = hdul[idx].header
    return header

def read_fits_data(fits_path: str, idx: int = 0) -> fits.PrimaryHDU:
    """
    Read the data of a FITS file.
    
    :param fits_path: Path to the FITS file
    :param idx: Index of the HDU to read the data from (default is 0 for the primary HDU)
    :return: FITS data as a numpy array
    """
    with fits.open(fits_path) as hdul:
        data = hdul[idx].data
    return data

def read_fits_data_shape(fits_path: str, idx: int = 0) -> tuple:
    """
    Read the shape of the data in a FITS file without loading the entire data into memory.
    
    :param fits_path: Path to the FITS file
    :param idx: Index of the HDU to read the data shape from (default is 0 for the primary HDU)
    :return: Shape of the FITS data
    """
    with fits.open(fits_path) as hdul:
        data_shape = hdul[idx].data.shape
    return data_shape
