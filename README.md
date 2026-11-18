# Skala-1.1 barriers and paths of five small molecules: code, data and companion website

This repository holds everything needed to look at, check and re-run the calculations of the study of reaction barriers and reaction paths of five small molecules (HCN to HNC, NH3 inversion, H2CO to H2 + CO, and two SN2 reactions) with the learned Skala-1.1 density functional, compared with coupled-cluster references. It contains the notebooks, the Python modules, all result tables, the figures and the script that draws them, and a static website that shows the numbers. It does not contain the text of the paper.

- Website (GitHub Pages): **PAGES_URL_PLACEHOLDER**
- Repository: `https://github.com/REPO_PLACEHOLDER`
- Licences: code MIT (`LICENSE`); data, tables and figures CC BY 4.0 (`LICENSE-DATA.txt`)

The placeholders `REPO_PLACEHOLDER` (written `OWNER/REPOSITORY`), `PAGES_URL_PLACEHOLDER`, `HOLDER_PLACEHOLDER`, `DOI_PLACEHOLDER`, `CONTACT_PLACEHOLDER` and `AUTHORS_PLACEHOLDER` are replaced by the owner after the repository is created (files: `README.md`, `LICENSE`, `LICENSE-DATA.txt`, `CITATION.cff`, `js/config.js`).

## What you can do here

| I want to ... | Do this | Needs |
| --- | --- | --- |
| look at the numbers | open the website, or the CSV files in `data/csv/` | a browser |
| re-check the published numbers | `python code/recheck_numbers.py` (a few seconds; expected: 139 items, `check_report.csv: identical to the published file`) | Python 3, numpy |
| redraw the figures and tables | `python figures/code/build_figures.py`, then `python figures/code/compare_with_published.py` (expected: `all agree`) | `requirements.txt` |
| re-run the calculations | the Colab notebooks below | a Google account; Colab CPU and GPU (T4) runtimes |
| rebuild the website pages | `python tools/build_site.py` and `python tools/check_links.py` | Python 3, numpy |

The website must be opened through a web server, not by double-clicking `index.html`; locally: `python -m http.server 8000`, then <http://localhost:8000/>.

## Reproduce the calculations on Colab

Run the CPU notebook first, then the GPU notebook. The notebooks mount your Google Drive (one click on **Allow**), write every result to a CSV file in `MyDrive/skala-b1-v5`, and can be started again after a disconnect (finished calculations are skipped). Expected times on a free runtime: CPU about 28 minutes (up to 50 with the tail), GPU about 50 minutes on a T4. Details and the SCF rescue ladder are on the website's Code page.

| Notebook | What it does | Colab |
| --- | --- | --- |
| `B1v5_CPU_1h_v2` | coupled-cluster reference side (CCSD(T), MP2 basis correction) | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/REPO_PLACEHOLDER/blob/main/code/notebooks/B1v5_CPU_1h_v2.ipynb) |
| `B1v5_GPU_1h_v3` | functionals: Skala-1.1, wB97M-V, B3LYP-D3(BJ), r2SCAN, PBE | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/REPO_PLACEHOLDER/blob/main/code/notebooks/B1v5_GPU_1h_v3.ipynb) |
| `B1v5_GPU_diag_sn2` | diagnostic 1: Cl- + CH3Cl with Skala-1.1 | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/REPO_PLACEHOLDER/blob/main/code/notebooks/B1v5_GPU_diag_sn2.ipynb) |
| `B1v5_GPU_diag_sym` | diagnostic 2: symmetric stretch | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/REPO_PLACEHOLDER/blob/main/code/notebooks/B1v5_GPU_diag_sym.ipynb) |
| `B1v5_collect_zip` | packs the results folder into one zip | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/REPO_PLACEHOLDER/blob/main/code/notebooks/B1v5_collect_zip.ipynb) |

The badges work after `REPO_PLACEHOLDER` is replaced by `OWNER/REPOSITORY` (edit this file with the pencil icon; the website's own buttons read `js/config.js`). The notebooks are stored without outputs. Software versions: `requirements-colab-cpu.txt` and `requirements-colab-gpu.txt` (pinned where the run pinned them; unpinned items are marked). The Skala-1.1 weights are downloaded from Hugging Face by the GPU notebook and are not part of this repository.

## How to cite

Cite the paper (details to be added on publication) and this repository. A machine-readable entry is in `CITATION.cff` (GitHub shows a "Cite this repository" box once the file is valid). The website's About page shows a citation box that fills in from `js/config.js`.

GitHub Pages is a web host, not an archive: a repository can be renamed, changed or deleted. A permanent copy with a DOI (for example on Zenodo) is a separate step; `SHA256SUMS.txt` lets anyone check that a copy and the website hold the same files.

## Licences

- Code (`code/`, `tools/`, `figures/code/`, `js/`, `css/`, the HTML pages, the notebooks): MIT, file `LICENSE`.
- Data, tables, figures and geometries (`data/`, `figures/png|pdf`, `b1_geometries.zip`): CC BY 4.0, file `LICENSE-DATA.txt` (full legal text at <https://creativecommons.org/licenses/by/4.0/legalcode>).
- Third-party software and model weights that the notebooks download (PySCF, GPU4PySCF, Skala and its weights, ...) have their own licences and are not distributed here.

## Folder map

    index.html  results/  demos/  figures/  data/  code/  checks/  about/    the eight pages of the website (each a folder with index.html)
    css/  js/                   styles and scripts; js/config.js holds the values to fill in (repository, holder, contact, DOI)
    data/csv/study|results|tables|figures/   all CSV tables (79 files); the data dictionary is on the data page
    data/json/                  the same numbers as JSON, read by the pages
    data/records/               small JSON records of the run, phase05_hcn_grid.npz
    data/geometries/            37 xyz files;  b1_geometries.zip  all of them in one zip
    figures/png|pdf/           the figures (main text, supplementary, summary graphic); no SVG copies in this package
    figures/code/               build_figures.py and its helper modules; compare_with_published.py; README.md
    code/notebooks/             five Colab notebooks, outputs removed
    code/modules/               the Python modules the notebooks paste into their cells, builders and small tests
    code/geompack/              the geometry pack (JSON) with its SHA-256
    code/recheck_numbers.py     re-runs the automatic check on the published records
    tools/                      builder of the pages, link checker, data dictionary
    requirements*.txt           Python packages (local use; Colab CPU; Colab GPU)
    SHA256SUMS.txt              SHA-256 of every file in data/, code/, figures/ and of b1_geometries.zip
    LICENSE  LICENSE-DATA.txt  CITATION.cff  .nojekyll  .gitignore

## Copy for double-blind review

Make a separate repository with a neutral name under an account that does not identify the authors, upload the same files, and set `anonymous: true` in `js/config.js`. The pages then carry `<meta name="robots" content="noindex">`, and the contact line, copyright holder, DOI and repository are left out of the citation box and the About page. Leave the holder, contact and author placeholders in this copy and fill them in after acceptance.

## About the image files (state of this package)

The environment in which the figures were produced added an invisible provenance block to every PNG and SVG file written in it (a `caBX` chunk in the PNG files, a `<metadata>` element in the SVG files). The block is not part of the pictures. Whether it is removed from the published files has not been decided. For now the files in the package folder are as written; the zip files of the package contain the files without the block (PNG chunks not needed to draw the image are dropped, the image data is copied byte for byte), and `SHA256SUMS.txt` lists the hashes of the files without the block, as in the earlier site package. If you upload a folder instead of the zip files and want `sha256sum -c SHA256SUMS.txt` to pass for the figures, first run `python tools/strip_png_metadata.py figures/png/*.png figures/svg/*.svg`. The pictures are identical either way.

## Known limits of the checks

- GitHub Pages itself was not tested from the build environment (github.io is not reachable from it); the pages were tested on a local web server with Chromium and with a link checker.
- `code/modules/test_1h_notebooks.py` needs PySCF and related packages and was not run for this package. `code/modules/check_test_synthetic.py` passes 100 assertions and fails one (it fails the same way in the unmodified source); `code/modules/test_1h_check.py` passes 57 of 57.
- The figure script gives pixel-identical PNG files with matplotlib 3.11.2 and Pillow 12.3.0; with other versions anti-aliasing can differ slightly.
