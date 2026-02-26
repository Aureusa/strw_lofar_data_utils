![Example LoTSS cutout](example_cutout.png)

# 🌌 STRW LOFAR Data Utils

This utility package was created out of necessity when no straightforward tools existed to generate cutouts from LoTSS (LOFAR Two-metre Sky Survey) DR2 data. It leverages the fact that most LoTSS DR2 mosaics are stored on the STRW cluster to provide fast and efficient cutout generation for any position within DR2 coverage.

Built primarily for STRW/Leiden Observatory environments.

## ✨ Features

- 📍 Generate cutouts from `(RA, Dec)` coordinates
- 📐 Specify size in pixels or arcminutes
- ⚙️ Batch processing via a single pipeline entrypoint
- 🧭 DR2 crawler support - tile-based cutout generation
- 📚 Source catalog integration - identify which sources fall within cutouts and get their positions. Source catalog here refers to `PyBDSF`-style catalogs containing source/component positions and properties.
- 🖼️ Cutout plotting and save-to-FITS helpers

## ✅ Requirements

- Python 3.9+
- Access to LoTSS DR2 mosaic files (cluster/local mirror)

## 📦 Installation

### Option 1 (recommended for development): editable install

```bash
git clone https://github.com/Aureusa/strw_lofar_data_utils.git
cd strw_lofar_data_utils

python -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e .
```

With editable install, code changes are picked up immediately (no reinstall needed for normal `.py` edits).

### Option 2: standard local install

```bash
python -m pip install .
```

## 🧩 Import Paths (current)

Use the real package path:

```python
from strw_lofar_data_utils.pipelines import generate_cutouts
from strw_lofar_data_utils.core.mosaic_crawler import DR2Crawler
from strw_lofar_data_utils.core.cutout_maker.cutout_catalogue import CutoutCatalogue
```

## 🚀 Usage
See the `examples/` directory for Jupyter notebooks demonstrating cutout generation, catalog integration, and visualization. All the main utilities in this package have been showcased in the notebooks with example code snippets and explanations.

## 🛠️ Environment Configuration (for STRW users you can ignore this as the shipped files already contains the correct paths/settings for STRW use)

The code uses `.env` values for cluster-specific paths/settings. Typical values:

```bash
BASE_DIR="/disks/paradata/shimwell/LoTSS-DR2/mosaics"
RA13h_field="RA13h_field"
RA0h_field="RA0h_field"
RMS="0.1"
```

If you run outside STRW, set `BASE_DIR` and related values to your local data layout. `Highly discuraged to run outside STRW` as the code assumptions are highly coupled to the speicific STRW environment. This might be resolved in the future by adding more configuration options and decoupling from STRW-specific paths but for now this package is developed `exclusively for STRW use`.

## 🗺️ First-Time Coverage File Generation (Not needed for STRW users as this is already done and the file is shipped with the package)

A precomputed coverage file is already shipped at:

- `data/mosaic_coverage/lotss_dr2_mosaic_coverage.csv`

Regenerate only if needed:

```bash
cd scripts
./get_mosaic_coverage.sh
```

## 🗂️ Project Layout

```text
strw_lofar_data_utils/
├── pyproject.toml
├── README.md
├── requirements.txt
├── data/ # precomputed files, e.g. mosaic coverage, catalogs
├── examples/ # Jupyter notebooks demonstrating usage
├── scripts/
└── strw_lofar_data_utils/
    ├── core/
    │   ├── cutout_maker/
    │   └── mosaic_crawler/
    └── pipelines/
```

## 📄 License

MIT. See [LICENSE](LICENSE).
