from astropy.io.fits.header import Header
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import pandas as pd

from ..utils import pixel_2_arcmin_size, arcmin_2_pixel_size


class AstroObject:
    def __init__(self, ra: list[float], dec: list[float]):
        self.ra = [r * u.deg for r in ra]
        self.dec = [d * u.deg for d in dec]

        self.x_pos = None
        self.y_pos = None

        self._in_cutout = None

        self._registered = False

    def get_coords(self) -> list[SkyCoord]:
        return [SkyCoord(ra=r, dec=d, frame='icrs') for r, d in zip(self.ra, self.dec)]

    def check_if_in_cutout(self, cutout_wcs, cutout_shape) -> list[bool]:
        if self._registered:
            return self._in_cutout
        in_cutout = []
        for coord in self.get_coords():
            pixel_coords = cutout_wcs.world_to_pixel(coord)
            x, y = int(pixel_coords[0]), int(pixel_coords[1])
            if 0 <= x < cutout_shape[1] and 0 <= y < cutout_shape[0]:
                in_cutout.append(True)
                self._set_pixel_positions(x, y)
            else:
                in_cutout.append(False)
        self._registered = True
        self._in_cutout = in_cutout
        return in_cutout
    
    def _set_pixel_positions(self, x: int, y: int) -> None:
        if self.x_pos is None:
            self.x_pos = []
        if self.y_pos is None:
            self.y_pos = []
        self.x_pos.append(x)
        self.y_pos.append(y)
        
    def get_pixel_positions(self) -> list[tuple[int, int]]:
        return list(zip(self.x_pos, self.y_pos))


def find_nearby_objects(
    ra_dec: tuple[float, float],
    header: Header,
    catalogue_path: str,
    size_pixels: int = None,
    size_arcmin: float = None,
    ) -> dict[str, AstroObject]:
    if size_pixels is None and size_arcmin is None:
        raise ValueError("Either size_pixels or size_arcmin must be provided.")
    
    if size_pixels is not None:
        max_sep_arcsec = pixel_2_arcmin_size(size_pixels, header) * 60 / 2
    else:
        max_sep_arcsec = size_arcmin * 60 / 2

    # Create a box around the given RA, Dec
    center_coord = SkyCoord(ra=ra_dec[0] * u.deg, dec=ra_dec[1] * u.deg, frame='icrs')
    delta_ra = max_sep_arcsec / np.cos(np.deg2rad(ra_dec[1])) / 3600  # in degrees
    delta_dec = max_sep_arcsec / 3600  # in degrees

    # Extract maximum and minimum RA and Dec
    ra_min = center_coord.ra.deg - delta_ra
    ra_max = center_coord.ra.deg + delta_ra
    dec_min = center_coord.dec.deg - delta_dec
    dec_max = center_coord.dec.deg + delta_dec

    # Here we assume tha catalogue is a CSV with 'RAJ2000' and 'DEJ2000' columns
    catalogue = pd.read_csv(catalogue_path)
    nearby_objects = catalogue[
        (catalogue['Component_RA'] >= ra_min) & (catalogue['Component_RA'] <= ra_max) &
        (catalogue['Component_DEC'] >= dec_min) & (catalogue['Component_DEC'] <= dec_max)
    ]

    # Here we assume that the 'recno' column uniquely identifies each object and group by it
    astro_objects = {}
    for _, group in nearby_objects.groupby('recno'):
        obj = AstroObject(ra=group['Component_RA'].tolist(), dec=group['Component_DEC'].tolist())
        astro_objects[group['recno'].iloc[0]] = obj
    return astro_objects
