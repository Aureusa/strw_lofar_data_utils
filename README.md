# STRW LOFAR Data Utils

A convenient toolkit for generating FITS cutouts from LoTSS DR2 mosaics at any RA/Dec position with customizable sizes. Built specifically for use on the STRW cluster at Leiden Observatory.

## Overview

This repository was created out of necessity when no straightforward tools existed to generate cutouts from LoTSS (LOFAR Two-metre Sky Survey) DR2 data. It leverages the fact that most LoTSS DR2 mosaics are stored on the STRW cluster to provide fast and efficient cutout generation for any position within DR2 coverage.

### Key Features

- **Simple Interface**: Provide RA/Dec coordinates and desired size - get your cutout
- **Flexible Sizing**: Specify cutout size in either arcminutes or pixels
- **Batch Processing**: Generate multiple cutouts efficiently with progress tracking
- **Visualization Tools**: Built-in plotting with customizable contour levels
- **Mosaic Discovery**: Automatic identification of the correct mosaic for any sky position
- **Coverage Mapping**: Pre-computed mosaic coverage for fast lookups

## Requirements

- Python 3.7+
- Access to STRW cluster (required for accessing LoTSS DR2 mosaics)
- Dependencies listed in `requirements.txt`

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Aureusa/strw_lofar_data_utils.git
cd strw_lofar_data_utils
```

2. Create and activate a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Linux/Mac
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

4. Configure environment variables:

The `.env` file contains paths and settings specific to the STRW cluster (generally these should not need to be changed):
```bash
# Base directory for LOFAR DR2 mosaics on the STRW cluster
BASE_DIR="/disks/paradata/shimwell/LoTSS-DR2/mosaics"

# List of RA fields in DR2
RA13h_field="RA13h_field"
RA0h_field="RA0h_field"

# RMS noise level in mJy/beam (used for contour levels)
RMS="0.1"
```

## Quick Start

### Generate Coverage Map (First Time Setup)

There is a csv file containing all the mosaic coverage information in `data/mosaic_coverage/lotss_dr2_mosaic_coverage.csv`. You can generate such a file by running the following script (again this shouldn't be necessary):

```bash
cd scripts
./get_mosaic_coverage.sh
```

### Basic Usage

```python
from src.pipelines import generate_cutouts

# Define positions (RA, Dec in degrees)
positions = [
    (180.5, 45.3),  # (RA, Dec)
    (185.2, 47.8),
    (190.1, 50.2)
]

# Generate cutouts with 300 pixel size
cutouts = generate_cutouts(
    ra_dec_list=positions,
    size_pixels=300,
    save=False  # Set to True to save FITS files
)

# Display a cutout with contours
cutouts[0].show_cutout(
    title="My Source",
    contour_levels=[3, 5, 10]  # Multiples of RMS
)
```

### Using Arcminute Size

```python
cutouts = generate_cutouts(
    ra_dec_list=positions,
    size_arcmin=10.0,  # 10 arcminute cutouts
    save=True,
    data_folder="./my_cutouts"
)
```

### Working with Catalogs

```python
import pandas as pd
from src.pipelines import generate_cutouts

# Load your catalog
catalog = pd.read_csv("my_sources.csv")

# Extract positions
positions = list(zip(catalog["RA"].values, catalog["Dec"].values))

# Generate cutouts
cutouts = generate_cutouts(
    ra_dec_list=positions,
    size_pixels=500,
    save=True,
    data_folder="./catalog_cutouts"
)

print(f"Generated {len(cutouts)} cutouts successfully!")
```

This will create a folder `catalog_cutouts` with FITS files for each source. The usual naming convention will be applied:
`<field_name>_cutout_RA<RA_str>_Dec<dec_str>_Size<cutout_size>px.fits`

## Project Structure

```
strw_lofar_data_utils/
├── .env                          # Environment configuration
├── requirements.txt              # Python dependencies
├── LICENSE                       # MIT License
├── README.md                     # This file
├── data/
│   └── mosaic_coverage/         # Mosaic coverage data
│       └── lotss_dr2_mosaic_coverage.csv
├── scripts/
│   └── get_mosaic_coverage.sh   # Generate mosaic coverage map
├── src/
│   ├── core/
│   │   ├── get_mosaic_coverage.py   # Coverage extraction script
│   │   ├── mosaic.py                # Mosaic class
│   │   └── cutout_maker/
│   │       ├── cutout.py            # Cutout class
│   │       ├── finder.py            # Mosaic finder utilities
│   │       └── make_cutout.py       # Cutout generation
│   └── pipelines/
│       └── make_cutouts.py          # High-level pipeline
└── examples/
    └── cutout_generation.ipynb      # Example notebook
```

## API Reference

### `generate_cutouts()`

Main function for generating cutouts.

**Parameters:**
- `ra_dec_list` (list): List of (RA, Dec) tuples in degrees
- `size_arcmin` (float, optional): Cutout size in arcminutes
- `size_pixels` (int, optional): Cutout size in pixels (provide either this or `size_arcmin`)
- `data_folder` (str, optional): Folder to save cutouts if `save=True`
- `save` (bool): Whether to save cutouts to FITS files (default: False)
- `mosaic_coverage_file` (str): Path to coverage file (default: 'default' uses bundled file)

**Returns:**
- List of `Cutout` objects

### `Cutout` Class

**Methods:**
- `show_cutout(title, contour_levels=None, cmap='inferno', colorbar=False, save=False, save_path=None)`: Display the cutout with optional contours
- `save_cutout(output_path, filename=None)`: Save cutout to FITS file
- `offload_data()`: Free memory by offloading data

**Properties:**
- `cutout`: Astropy `Cutout2D` object
- `ra`, `dec`: Position coordinates
- `size_arcmin`, `size_pixels`: Cutout dimensions
- `mosaic`: Parent `Mosaic` object

### `Mosaic` Class

**Methods:**
- `is_in_coverage(ra, dec)`: Check if coordinates are within mosaic coverage
- `load_data()`: Load FITS data into memory
- `load_header()`: Load FITS header
- `offload_data()`: Free memory

**Properties:**
- `header`: FITS header
- `data`: FITS data array (must call `load_data()` first)
- `ra`, `dec`: Mosaic center coordinates
- `ra_min`, `ra_max`, `dec_min`, `dec_max`: Coverage boundaries

## Examples

See the `examples/cutout_generation.ipynb` notebook for a complete walkthrough including:
- Loading source catalogs
- Batch cutout generation
- Visualization with custom contours
- Memory management for large batches

## Configuration

### Mosaic Coverage Script

Edit `scripts/get_mosaic_coverage.sh` to configure:
- `PROJECTPATH`: Path to this repository
- `OUTPUTPATH`: Where to save the coverage CSV
- `SAVE`: Whether to save output (true/false)
- `VERBOSE`: Print detailed progress (true/false)

### Environment Variables

The `.env` file defines:
- `BASE_DIR`: Root directory of LoTSS DR2 mosaics on cluster
- `RA0h_field`, `RA13h_field`: Field naming conventions
- `RMS`: Default RMS noise level for contours (mJy/beam)

## Limitations

- **STRW Cluster Only**: This tool assumes specific file paths and naming conventions used on the STRW cluster at Leiden Observatory
- **DR2 Coverage**: Only works for positions within LoTSS DR2 coverage
- **Field Structure**: Assumes standard LoTSS DR2 field organization (RA0h_field, RA13h_field)

## Contributing

Pull requests and ideas to improve this tool are welcome! If you encounter issues or have suggestions:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-idea`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/your-idea`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Built for internal use at Leiden Observatory (STRW)
- Uses data from the LOFAR Two-metre Sky Survey (LoTSS) DR2
- Thanks to the LoTSS team for making the data accessible on the cluster

## Contact

For questions or issues specific to the STRW cluster setup, please contact the repository maintainer.

---

**Note**: This tool is designed for internal use at STRW and requires access to the Leiden Observatory computing cluster.
