# Figure build

`code/build_figures.py` rebuilds every figure, table, caption and check of the paper from the files in `data/` of this repository.

    pip install -r requirements.txt                    # from the repository root; numpy, pandas, matplotlib, pillow
    python figures/code/build_figures.py               # 20-30 seconds; writes figures_rebuilt/
    python figures/code/compare_with_published.py      # compares with data/csv/figures, data/csv/tables, figures/png

Paths are found from the location of the script, so the working directory does not matter. The published figures in `figures/png|pdf|svg` are never overwritten; the rebuild goes to `figures_rebuilt/` (ignored by git), or to a folder named with `--out`.

| module | content |
| --- | --- |
| `build_figures.py` | entry point: runs everything, writes `CAPTIONS.md`, `FIGURE_BUILD_LOG.md` and the colour-vision checks |
| `figcommon.py` | style, method colours, paths, data loading |
| `figs_fallback.py`, `figs_rest.py`, `figs_extra.py` | the figures (main text, then supplementary) |
| `tables.py` | Tables 1 to 5 and S1 to S8 as CSV and Markdown |
| `checks.py` | grey-scale and deuteranopia copies, contact sheets |
| `captions.py` | figure captions |
| `visual_qa_notes.md` | notes from the visual check, copied into the build log |

Input files (all inside this repository): `data/csv/study/` (master_*, summary_stats, path_profiles, sn2_symscan, derived_*), `data/csv/results/` (phase05_hcn_grid.csv, phase11_dft.csv, skala_d3_table.csv, table1_systems.csv, findings_outcome.csv, geometry_table.csv, excluded_and_failed.csv), `data/records/` (phase05_hcn_grid.npz, phase13_xtb_barriers.json), `data/geometries/` (GFN2-xTB transition-state xyz files).

Text strings in the table captions and in the build log that read `W/...` mean `data/csv/study/...`, and `R/...` mean `data/csv/results/...` (or `data/records/`); they are labels from the original run and are left unchanged so that the rebuilt files are identical to the published ones.

Expected result: 80 anchor checks, 79 pass and 1 is a known, documented mismatch (Figure 8: the plan wording "only the H2CO transition state differs by more than 0.35" gives 0.43 for NH3; the caption uses 0.5 instead, see `FIGURE_BUILD_LOG.md`). The PNG files are pixel-identical to the published ones with matplotlib 3.11.2 and Pillow 12.3; other versions may differ by a few pixels in anti-aliasing and font hinting, the CSV files do not depend on them.
