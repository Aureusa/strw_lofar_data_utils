from abc import ABC, abstractmethod


class BaseMosaicCrawler(ABC):
    def __post_init__(self):
        if not hasattr(self, '_x_list'):
            error_msg = "Subclasses of BaseMosaicCrawler must have an '_x_list' attribute."
            error_msg += " This list should contain the x pixel coordinates of the cutouts generated during crawling."
            raise AttributeError(error_msg)
        if not hasattr(self, '_y_list'):
            error_msg = "Subclasses of BaseMosaicCrawler must have a '_y_list' attribute."
            error_msg += " This list should contain the y pixel coordinates of the cutouts generated during crawling."
            raise AttributeError(error_msg)
        if not hasattr(self, '_cutout_size'):
            error_msg = "Subclasses of BaseMosaicCrawler must have a '_cutout_size' attribute."
            error_msg += " This should be an integer representing the size of the cutouts generate duiring the crawling in pixels."
            raise AttributeError(error_msg)
        if not hasattr(self, '_mosaic_shape'):
            error_msg = "Subclasses of BaseMosaicCrawler must have a '_mosaic_shape' attribute."
            error_msg += " This should be a tuple representing the shape of the mosaic in pixels (width, height)."
            raise AttributeError(error_msg)

    @abstractmethod
    def crawl(self):
        """
        Abstract method to crawl through the mosaic and generate cutouts based on the RA and Dec list.
        """
        pass

    @abstractmethod
    def has_crawled(self) -> bool:
        """
        Abstract method to check if the crawling process has been completed.
        """
        pass

    @property
    def x_list(self) -> list:
        """
        Get the list of x pixel coordinates of the cutouts generated during crawling.

        :return: List of x pixel coordinates.
        """
        return self._x_list
    
    @property
    def y_list(self) -> list:
        """
        Get the list of y pixel coordinates of the cutouts generated during crawling.

        :return: List of y pixel coordinates.
        """
        return self._y_list
    
    @property
    def cutout_size(self) -> int:
        """
        Get the size of the cutouts in pixels.

        :return: Size of the cutouts in pixels.
        """
        return self._cutout_size
    
    @property
    def mosaic_shape(self) -> tuple:
        """
        Get the shape of the mosaic in pixels.

        :return: Shape of the mosaic in pixels (width, height).
        """
        return self._mosaic_shape