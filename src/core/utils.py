from astropy.io.fits.header import Header


def pixel_2_arcmin_size(size_pixels: int, header: Header) -> float:
        """
        Calculate size in arcminutes from size in pixels.
        
        :param size_pixels: Size in pixels
        :return: Size in arcminutes
        """
        pixel_scale_deg = abs(header['CDELT1'])  # degrees per pixel
        size_deg = size_pixels * pixel_scale_deg
        return size_deg * 60  # convert to arcminutes
    
def arcmin_2_pixel_size(size_arcmin: float, header: Header) -> int:
    """
    Calculate size in pixels from size in arcminutes.
    
    :param size_arcmin: Size in arcminutes
    :return: Size in pixels
    """
    pixel_scale_deg = abs(header['CDELT1'])  # degrees per pixel
    size_deg = size_arcmin / 60  # convert to degrees
    return int(size_deg / pixel_scale_deg)
