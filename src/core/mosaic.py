import os
import glob
import dotenv
import gc

from astropy.io import fits
import pandas as pd


dotenv.load_dotenv()
BASE_DIR = os.getenv("BASE_DIR", "/disks/paradata/shimwell/LoTSS-DR2/mosaics")
RA0h_field = os.getenv("RA0h_field", "RA0h_field")
RA13h_field = os.getenv("RA13h_field", "RA13h_field")
FIELD_LIST = [RA0h_field, RA13h_field]


class Mosaic:
    def __init__(
            self,
            field_name: str,
            ra: float, dec: float,
            ra_min: float, 
            ra_max: float,
            dec_min: float,
            dec_max: float,
            ra_size: float,
            dec_size: float
        ):
        self.field_name = field_name
        self.ra = ra
        self.dec = dec
        self.ra_min = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max
        self.ra_size = ra_size
        self.dec_size = dec_size

        self._header = None
        self._data = None

        # Find the mosaic paths
        self._field_path = self._find_field_path()
        self._mosaic_blanked_path = self._find_mosaic_blanked_path()

    @property
    def header(self):
        if self._header is None:
            print("Header not loaded. Call load_header() to load it first.")
        return self._header
    
    @property
    def data(self):
        if self._data is None:
            print("Data not loaded. Call load_data() to load it first.")
        return self._data

    def is_in_coverage(self, ra: float, dec: float) -> bool:
        """
        Check if the given RA and Dec are within the mosaic coverage.

        :param ra: Right Ascension in degrees
        :param dec: Declination in degrees
        :return: True if within coverage, False otherwise
        """
        return (self.ra_min <= ra <= self.ra_max) and (self.dec_min <= dec <= self.dec_max)
    
    def load_header(self):
        """
        Load the FITS header for the mosaic.

        :return: FITS header
        """
        if self._header is None:
            fits_path = self._mosaic_blanked_path
            with fits.open(fits_path) as hdul:
                self._header = hdul[0].header
        return self._header
    
    def load_data(self):
        """
        Load the FITS data for the mosaic.

        :return: FITS data
        """
        if self._data is None:
            fits_path = self._mosaic_blanked_path
            with fits.open(fits_path) as hdul:
                self._data = hdul[0].data
        return self._data
    
    def offload_data(self):
        """
        Offload the FITS data from memory.
        """
        self._data = None
        self._header = None
        gc.collect()

    def _find_mosaic_blanked_path(self):
        mosaic_blanked = "mosaic-blanked.fits"
        return os.path.join(
            self._field_path,
            self.field_name,
            mosaic_blanked
        )

    def _find_field_path(self) -> str:
        # Get all subdirectories in BASE_DIR + the two field names
        search_paths = [
            os.path.join(BASE_DIR, RA0h_field),
            os.path.join(BASE_DIR, RA13h_field)
        ]

        for base_path in search_paths:
            subdirs = glob.glob(os.path.join(base_path, "*"))
            for subdir in subdirs:
                if os.path.basename(subdir) == self.field_name:
                    return base_path
        raise FileNotFoundError(f"Field {self.field_name} not found in expected directories.")
    

def get_list_of_mosaics(mosaic_coverage_file: str = 'default') -> list[Mosaic]:
    if mosaic_coverage_file == 'default':
        mosaic_coverage_file = os.path.join(
            os.path.dirname(__file__),
            '..',
            '..',
            'data',
            'mosaic_coverage',
            'lotss_dr2_mosaic_coverage.csv'
        )
    df = pd.read_csv(mosaic_coverage_file)

    mosaics = []
    for _, row in df.iterrows():
        mosaic = Mosaic(
            field_name=row['field_name'],
            ra=row['ra'],
            dec=row['dec'],
            ra_min=row['ra_min'],
            ra_max=row['ra_max'],
            dec_min=row['dec_min'],
            dec_max=row['dec_max'],
            ra_size=row['ra_size'],
            dec_size=row['dec_size']
        )
        mosaics.append(mosaic)
    return mosaics
    