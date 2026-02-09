from astropy.visualization import (AsinhStretch, ImageNormalize, PercentileInterval)
import numpy as np
import matplotlib.pyplot as plt

from .cutout_maker.cutout import Cutout
from .mosaic_crawler.base_crawler import BaseMosaicCrawler


class Vizualizer:
    def vizualize_a_cutout(self, cutout: Cutout, title: str, contour_levels: list = None, cmap='inferno', colorbar=False, rms=0.1, save_path=None, save=False) -> None:
        """
        Show the cutout using matplotlib.

        :param title: Title for the plot
        :param contour_levels: List of contour levels in multiples of RMS noise (e.g., [3, 5, 10])
        :param cmap: Colormap for the image
        :param colorbar: Whether to show a colorbar
        :param rms: RMS noise level in mJy/beam
        (used to calculate contour levels if contour_levels is provided)
        :param save_path: Path to save the plot image (if save is True)
        :param save: Whether to save the plot instead of showing it
        """
        cutout_wcs = cutout.get_wcs()
        data = cutout.get_data()

        contour_levels = self._make_contour_levels(contour_levels, rms) if contour_levels is not None else None

        norm = ImageNormalize(data, interval=PercentileInterval(99.5), stretch=AsinhStretch())
        fig = plt.figure(figsize=(10, 10))
        ax = plt.subplot(projection=cutout_wcs)
        im = ax.imshow(data, cmap=cmap, origin='lower', norm=norm)
        if colorbar:
            plt.colorbar(im, ax=ax, label='Intensity')
        if contour_levels is not None:
            ax.contour(data, levels=contour_levels, colors='white', linewidths=0.5)
        ax.set_title(title)
        ax.set_xlabel('RA (J2000)')
        ax.set_ylabel('Dec (J2000)')
        if save:
            if save_path is None:
                save_path = f"{title.replace(' ', '_')}_cutout.png"
            plt.savefig(save_path)
        else:
            plt.show()
        plt.close()

    def vizualize_a_crawl(self, crawler: BaseMosaicCrawler) -> None:
        if not crawler.has_crawled():
            raise ValueError("The crawler has not completed crawling yet. Please call the 'crawl' method before visualizing.")
        mosaic_radius = crawler.mosaic_shape[0] / 2
        crawler_cutout_size = crawler.cutout_size
        mosaic_wcs = crawler.mosaic.wcs

        fig = plt.figure(figsize=(8, 8))
        ax = plt.subplot(projection=mosaic_wcs)
        ax.set_title(f"Crawl Visualization of Mosaic: {crawler.mosaic.field_name}")
        ax.set_xlabel('RA (J2000)')
        ax.set_ylabel('Dec (J2000)')

        # Plot the mosaic coverage as a circle
        circle = plt.Circle((mosaic_radius, mosaic_radius), mosaic_radius, color='lightgray', alpha=0.5)
        ax.add_artist(circle)

        # In the top left corner add a rectangle representing the cutout size
        rect = plt.Rectangle((10, 10), crawler_cutout_size, crawler_cutout_size, edgecolor='red', facecolor='none', alpha=0.5)
        ax.add_artist(rect)
        ax.text(10, 10 + crawler_cutout_size + 20, f"Cutout Size: {crawler_cutout_size} pixels", color='red')

        # Scatter the positions of the generated cutouts
        x_array = np.array(crawler.x_list)
        y_array = np.array(crawler.y_list)
        ax.scatter(x_array, y_array, color='blue', s=10, label='Generated Cutouts Centers')

        # Plot rectangles representing the cutout areas
        for x, y in zip(x_array, y_array):
            rect = plt.Rectangle(
                (x - crawler_cutout_size / 2, y - crawler_cutout_size / 2),
                crawler_cutout_size,
                crawler_cutout_size,
                edgecolor='red',
                facecolor='none',
                alpha=0.5
            )
            ax.add_artist(rect)
        ax.legend()
        plt.xlim(0, crawler.mosaic_shape[0])
        plt.ylim(0, crawler.mosaic_shape[1])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def _make_contour_levels(self, levels: list, rms: float) -> np.ndarray:
        """
        Generates contour levels based on RMS noise.
        
        :param levels: List of contour levels in multiples of RMS noise
        :return: List of contour levels in Jy/beam
        """
        return np.array(levels) * rms * 1e-3  # Convert mJy/beam to Jy/beam