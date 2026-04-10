import numpy as np
from tqdm import tqdm

from .base_crawler import BaseMosaicCrawler
from ..mosaic import get_mosaic_by_field_name
from ..cutout_maker import make_cutout, Cutout
from ...logger import setup_logging

class Crawler(BaseMosaicCrawler):
    def __init__(
            self,
            field_name: str,
            cutout_size: int,
            stride: float,
            mosaic_coverage_file: str = 'default',
            release: str = 'DR2',
            verbose: bool = False
        ):
        """
        Initialize the Crawler with the given parameters.
        
        :param field_name: Name of the mosaic field to crawl
        :type field_name: str
        :param cutout_size: Size of the cutouts to generate in pixels
        :type cutout_size: int
        :param stride: Stride in decimal percentage. For example:
            stride=0.5: the cutouts will be generated with a 50% overlap.
            stride=0.1: the cutouts will be generated with a 90% overlap.
            stride=1: the cutouts will be generated without overlap.
            stride=2: the cutouts will be generated with a gap of 100%
                of the cutout size between them.
        :type stride: float
        :param mosaic_coverage_file: Path to the mosaic coverage CSV file. If 'default'
        uses the coverage file in the `/path/to/repo/data/mosaic_coverage/lotss_dr2_mosaic_coverage.csv`.
        :type mosaic_coverage_file: str
        :param release: Data release to use ('DR2' or 'DR3')
        :type release: str
        """
        # Get the mosaic for the given field name and its shape
        self.mosaic = get_mosaic_by_field_name(
            field_name, mosaic_coverage_file=mosaic_coverage_file, release=release
        )
        self._mosaic_shape = self.mosaic.get_data_shape()

        self._cutout_size = cutout_size
        self.stride = stride

        num_cutouts_x = self._mosaic_shape[0] // (cutout_size * stride)
        num_cutouts_y = self._mosaic_shape[1] // (cutout_size * stride)
        self.total_cutouts = int(num_cutouts_x * num_cutouts_y)


        self.ra_dec_list, self._x_list, self._y_list = self._generate_ra_dec_list(
            num_cutouts_x, num_cutouts_y, self._cutout_size, self.stride
        )

        self.verbose = verbose

        if verbose:
            self.logger = setup_logging(name="strw_lofar_data_utils.core.mosaic_crawler.crawler.Crawler")
            self._log_info()

    def crawl(self) -> list[Cutout]:
        """
        Crawl through the mosaic and generate cutouts based on the RA and Dec list.
        """
        new_ra_dec_list = []
        new_x_list = []
        new_y_list = []

        cutouts = []
        for idx, (ra, dec) in enumerate(
            tqdm(self.ra_dec_list, desc="Generating cutouts") if self.verbose
            else self.ra_dec_list
            ):
            cutout = make_cutout(self.mosaic, ra, dec, size_pixels=self._cutout_size)
            if cutout is not None:
                cutouts.append(cutout)
                new_ra_dec_list.append((ra, dec))
                new_x_list.append(self._x_list[idx])
                new_y_list.append(self._y_list[idx])
        if self.verbose:
            self.logger.info(f"Generated {len(cutouts)} cutouts out of {len(self.ra_dec_list)} possible positions.")

        self.ra_dec_list = new_ra_dec_list
        self._x_list = new_x_list
        self._y_list = new_y_list

        self._has_crawled = True
        return cutouts
    
    def has_crawled(self) -> bool:
        """
        Check if the crawling process has been completed.

        :return: True if crawling is completed, False otherwise.
        """
        return hasattr(self, '_has_crawled') and self._has_crawled

    def _log_info(self):
        self.logger.info(f"Mosaic field name: {self.mosaic.field_name}")
        self.logger.info(f"Cutout size (pixels): {self._cutout_size}")
        self.logger.info(f"Stride: {self.stride}")
        self.logger.info(f"Total cutouts that can be generated assuming square mosaic ({self._mosaic_shape}): {self.total_cutouts}")
        self.logger.info(f"Warning: This is an upper limit on the number of cutouts that can be generated.")
        self.logger.info("The actual number will be lower due to edge effects and the shape of the mosaic coverage (circle).")

    def _generate_ra_dec_list(self, num_cutouts_x, num_cutouts_y, cutout_size, stride):
        # Generate indices using numpy for vectorized operation
        i_indices = np.arange(num_cutouts_x)
        j_indices = np.arange(num_cutouts_y)
        
        # Create meshgrid for all combinations of i and j
        i_grid, j_grid = np.meshgrid(i_indices, j_indices, indexing='ij')
        
        # Calculate x and y positions vectorized
        x_arr = (i_grid * cutout_size * stride).flatten()
        y_arr = (j_grid * cutout_size * stride).flatten()

        x_list = x_arr.tolist()
        y_list = y_arr.tolist()

        ra_dec_list = self.mosaic.convert_pixels_to_world(x_list, y_list)
        return ra_dec_list, x_list, y_list
    