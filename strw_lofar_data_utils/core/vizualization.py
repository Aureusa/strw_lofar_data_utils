from astropy.visualization import (
    AsinhStretch,
    ImageNormalize,
    PercentileInterval,
    BaseStretch,
    BaseInterval
)
import numpy as np
import matplotlib.pyplot as plt

from .cutout_maker import Cutout, CutoutCatalogue
from .mosaic_crawler.base_crawler import BaseMosaicCrawler
from .mosaic import Mosaic


class Vizualizer:
    def vizualize_a_cutout(
            self,
            cutout: Cutout,
            title: str,
            contour_levels: list = None,
            interval: BaseInterval = PercentileInterval(99.5),
            stretch: BaseStretch = AsinhStretch(),
            cmap='inferno',
            colorbar=False,
            rms=0.1,
            save_path=None,
            save=False
        ) -> None:
        """
        Show the cutout using matplotlib.

        :param cutout: Cutout object to visualize
        :param title: Title for the plot
        :param contour_levels: List of contour levels in multiples of RMS noise (e.g., [3, 5, 10])
        :param interval: Interval object for image normalization
        :param stretch: Stretch object for image normalization
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

        norm = self._get_norm(data, interval=interval, stretch=stretch)
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

    def vizualize_a_cutout_with_catalogue_objects(
            self,
            cutout: Cutout,
            cutout_cat: CutoutCatalogue,
            title: str,
            contour_levels: list = [4, 6, 8, 10],
            interval: BaseInterval = PercentileInterval(99.5),
            stretch: BaseStretch = AsinhStretch(),
            cmap='gray',
            rms=0.1,
            save_path=None,
            legend=True,
            save=False
        ) -> None:
        """
        Visualize a cutout with overlaid positions of objects from the cutout catalogue.
        
        :param cutout: Cutout object to visualize
        :param cutout_cat: CutoutCatalogue object containing the catalog of objects in the cutout
        :param title: Title for the plot
        :param contour_levels: List of contour levels in multiples of RMS noise (e.g., [3, 5, 10])
        :param interval: Interval object for image normalization
        :param stretch: Stretch object for image normalization
        :param cmap: Colormap for the image
        :param rms: RMS noise level in mJy/beam (used to calculate contour levels if contour_levels is provided)
        :param save_path: Path to save the plot image (if save is True)
        :param legend: Whether to show a legend for the catalog objects
        :param save: Whether to save the plot instead of showing it
        """
        data_cutout = cutout.get_data()

        # Creates a dict of SourceBlob instances for each unique object in the cutout
        sb_dict = cutout_cat.get_source_blobs_from_catalogue()

        # For plotting purposes
        norm = self._get_norm(data_cutout, interval=interval, stretch=stretch)
        if contour_levels is not None:
            contours = self._make_contour_levels(contour_levels, rms=rms)
        else:
            contours = None

        # Show the cutout
        _, ax = plt.subplots(figsize=(10, 10))
        ax.imshow(data_cutout, origin='lower', cmap=cmap, norm=norm)
        ax.contour(
            data_cutout,
            levels=contours,
            colors='black',
            linewidths=0.5,
            alpha=0.7
        )

        # Plot the positions of the objects within the cutout
        # as circular markers with different colors for each object
        colors = plt.cm.tab10(range(len(sb_dict)))
        for idx, (obj_id, source_blob) in enumerate(sb_dict.items()):
            sb_pixel_loc = source_blob.get_pixel_positions()
            x_coords = [pos[0] for pos in sb_pixel_loc] # Extract x coordinates for plotting
            y_coords = [pos[1] for pos in sb_pixel_loc] # Extract y coordinates for plotting
            ax.scatter(
                x_coords,
                y_coords,
                edgecolor=colors[idx],
                facecolor='none',
                s=5,
                linewidths=2,
                label=f"{idx} - {obj_id}"
            )
            for x, y in zip(x_coords, y_coords):
                ax.text(
                    x,
                    y,
                    str(idx),
                    color=colors[idx],
                    fontsize=24,
                )

        if legend:
            ax.legend(loc='upper right', fontsize='small')
        ax.set_title("Cutout with Catalog Objects Overlay")
        if save:
            if save_path is None:
                save_path = f"{title.replace(' ', '_')}_cutout_with_catalogue.png"
            plt.savefig(save_path)
        else:
            plt.show()
        plt.close()

    def vizualize_a_crawl(self, crawler: BaseMosaicCrawler) -> None:
        """
        Visualize the crawl process by plotting the positions of the generated cutouts on the mosaic.
        
        :param crawler: BaseMosaicCrawler object containing the crawl information
        """
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

    def vizualize_mosaic_coverage(self, mosaics: list[Mosaic]) -> None:
        """
        Visualize the coverage of multiple mosaics in RA/Dec space.
        Each mosaic is represented as a circle centered on its RA/Dec with a radius corresponding to its size.
        
        :param mosaics: List of Mosaic objects to visualize
        """
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='aitoff')
        ax.set_title("Mosaic Coverage in RA/Dec")
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')

        xticks_deg = np.arange(-150, 181, 30)
        ax.set_xticks(np.radians(xticks_deg))
        ax.set_xticklabels([f"{int((180 - deg) % 360)}°" for deg in xticks_deg])

        for mosaic in mosaics:
            ra_center_deg = mosaic.ra
            dec_center_deg = mosaic.dec

            # Keep original sky coordinates and transform only for Aitoff display,
            # centered at RA=180° and with astronomy-style RA direction.
            ra_center_shifted_deg = ((ra_center_deg - 180.0 + 360.0) % 360.0) - 180.0
            ra_center_rad = np.radians(-ra_center_shifted_deg)
            dec_center_rad = np.radians(dec_center_deg)
            ra_size_rad = np.radians(mosaic.ra_size)

            # Plot the center of the mosaic
            ax.scatter(ra_center_rad, dec_center_rad, marker='o', label=mosaic.field_name if len(mosaics) <= 10 else None)

            # Plot the coverage area as a circular patch (approximating the rectangular coverage in RA/Dec)
            circle = plt.Circle((ra_center_rad, dec_center_rad), ra_size_rad / 2, edgecolor='blue', facecolor='none', alpha=0.5)
            ax.add_artist(circle)

        if len(mosaics) <= 10:
            ax.legend(loc='upper right')
        plt.grid()
        plt.show()

    def _make_contour_levels(self, levels: list, rms: float) -> np.ndarray:
        """
        Generates contour levels based on RMS noise.
        
        :param levels: List of contour levels in multiples of RMS noise
        :param rms: RMS noise value in mJy/beam
        :return: List of contour levels in Jy/beam
        """
        return np.array(levels) * rms * 1e-3  # Convert mJy/beam to Jy/beam

    def _get_norm(self, data: np.ndarray, interval: BaseInterval, stretch: BaseStretch) -> ImageNormalize:
        """
        Get an ImageNormalize object for the given data using an interval and stretch.
        
        :param data: 2D numpy array of the image data
        :param interval: Interval object for image normalization
        :param stretch: Stretch object for image normalization
        :return: ImageNormalize object for normalizing the image display
        """
        return ImageNormalize(data, interval=interval, stretch=stretch)
    