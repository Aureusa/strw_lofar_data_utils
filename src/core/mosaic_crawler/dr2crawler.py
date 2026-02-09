import numpy as np
from tqdm import tqdm

from .base_crawler import BaseMosaicCrawler
from ..mosaic import get_mosaic_by_field_name
from ..cutout_maker import make_cutout


class DR2Crawler(BaseMosaicCrawler):
    def __init__(
            self,
            field_name: str,
            cutout_size: int,
            stride: float,
        ):
        """
        Initialize the DR2Crawler with the given parameters.
        
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
        """
        # Get the mosaic for the given field name and its shape
        self.mosaic = get_mosaic_by_field_name(field_name)
        self._mosaic_shape = self.mosaic.get_data_shape()

        self._cutout_size = cutout_size
        self.stride = stride

        num_cutouts_x = self._mosaic_shape[0] // (cutout_size * stride)
        num_cutouts_y = self._mosaic_shape[1] // (cutout_size * stride)
        self.total_cutouts = int(num_cutouts_x * num_cutouts_y)


        self.ra_dec_list, self._x_list, self._y_list = self._generate_ra_dec_list(
            num_cutouts_x, num_cutouts_y, self._cutout_size, self.stride
        )

        self._log_info()

    def crawl(self):
        """
        Crawl through the mosaic and generate cutouts based on the RA and Dec list.
        """
        new_ra_dec_list = []
        new_x_list = []
        new_y_list = []

        cutouts = []
        for idx, (ra, dec) in enumerate(tqdm(self.ra_dec_list, desc="Generating cutouts")):
            cutout = make_cutout(self.mosaic, ra, dec, size_pixels=self._cutout_size)
            if cutout is not None:
                cutouts.append(cutout)
                new_ra_dec_list.append((ra, dec))
                new_x_list.append(self._x_list[idx])
                new_y_list.append(self._y_list[idx])
        print(f"Generated {len(cutouts)} cutouts out of a possible {self.total_cutouts} based on the RA and Dec list.")

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
    
    def visualize_the_crawl(self):
        import matplotlib.pyplot as plt

        mosaic_radius = self._mosaic_shape[0] / 2

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_title(f"Crawl Visualization for Mosaic: {self.mosaic.field_name}")
        ax.set_xlabel('X Pixel')
        ax.set_ylabel('Y Pixel')

        # Plot the mosaic coverage as a circle
        circle = plt.Circle((mosaic_radius, mosaic_radius), mosaic_radius, color='lightgray', alpha=0.5)
        ax.add_artist(circle)

        # In the top left corner add a rectangle representing the cutout size
        rect = plt.Rectangle(10, 10, self._cutout_size, self._cutout_size, edgecolor='red', facecolor='none', alpha=0.5)
        ax.add_artist(rect)
        ax.text(10, 10 + self._cutout_size + 20, f"Cutout Size: {self._cutout_size} pixels", color='red')

        # Scatter the positions of the generated cutouts
        x_array = np.array(self._x_list)
        y_array = np.array(self._y_list)
        ax.scatter(x_array, y_array, color='blue', s=10, label='Generated Cutouts Centers')

        # Plot rectangles representing the cutout areas
        for x, y in zip(x_array, y_array):
            rect = plt.Rectangle(
                (x - self._cutout_size / 2, y - self._cutout_size / 2),
                self._cutout_size,
                self._cutout_size,
                edgecolor='red',
                facecolor='none',
                alpha=0.5
            )
            ax.add_artist(rect)
        ax.legend()
        plt.xlim(0, self._mosaic_shape[0])
        plt.ylim(0, self._mosaic_shape[1])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def _log_info(self):
        print(f"Mosaic field name: {self.mosaic.field_name}")
        print(f"Cutout size (pixels): {self._cutout_size}")
        print(f"Stride: {self.stride}")
        print(f"Total cutouts that can be generated assuming square mosaic ({self._mosaic_shape}): {self.total_cutouts}")
        print(f"Warning: This is an upper limit on the number of cutouts that can be generated.")
        print("The actual number will be lower due to edge effects and the shape of the mosaic coverage (circle).")

    def _generate_ra_dec_list(self, num_cutouts_x, num_cutouts_y, cutout_size, stride):
        # Generate indices using numpy for vectorized operation
        i_indices = np.arange(num_cutouts_x)
        j_indices = np.arange(num_cutouts_y)
        
        # Create meshgrid for all combinations of i and j
        i_grid, j_grid = np.meshgrid(i_indices, j_indices, indexing='ij')
        
        # Calculate x and y positions vectorized
        x_list = (i_grid * cutout_size * stride).flatten()
        y_list = (j_grid * cutout_size * stride).flatten()

        ra_dec_list = self.mosaic.convert_pixels_to_world(x_list, y_list)
        return ra_dec_list, x_list, y_list
    