from pathlib import Path

from astropy.io.fits.header import Header


def find_repo_root() -> Path:
        """
        Find the cloned repository root for this project.

        The package currently expects the top-level ``data/`` directory and ``.env``
        file to live in the repository checkout, not inside installed package data.
        """
        candidate_starts = [Path.cwd(), Path(__file__).resolve()]
        seen: set[Path] = set()

        for start in candidate_starts:
                for candidate in (start, *start.parents):
                        if candidate in seen:
                                continue
                        seen.add(candidate)

                        if (
                                (candidate / "pyproject.toml").exists()
                                and (candidate / "data").is_dir()
                                and (candidate / "strw_lofar_data_utils").is_dir()
                        ):
                                return candidate

        return None


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
