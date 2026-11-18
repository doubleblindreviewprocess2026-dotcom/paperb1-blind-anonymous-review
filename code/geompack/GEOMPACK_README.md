# Geometry pack of the 1-hour notebooks (b1v5-geompack-1)

File: `b1v5_geompack_1h.json`  (57 KB canonical JSON)  
SHA-256: `81e4533319f5fc593cc3e81daf93e6e4db4998df01eea275b3c813fc87dafae6`  
Built: 2026-11-15T14:56:39 in the sandbox; design: 1h core: six reactions + CH4 + ext78 (F- + CH3Cl, OH- + CH3F)  
Systems: hcn, nh3, h2co, ch3nc, sn2_ff, sn2_clcl, ch4, sn2_fcl, sn2_ohf

## What it contains

Every geometry the two 1-hour notebooks need, computed ONCE with the helper functions of the full-size notebooks (`b1v5_cpu.py`, driven by `b1v5_geompack.py`):

* optimised minima and transition states (own protocol: DF-B3LYP-D3(BJ), geomeTRIC, `GAU` convergence, tightened once with `GAU_TIGHT` if the RMS gradient exceeds 1.5e-4 hartree/bohr);
* harmonic frequencies from the analytic Hessian, the number of imaginary modes (transition states: one wanted; if more, one retry along the extra mode, then the scan peak is stored as a flagged fallback), zero-point energies;
* the relaxed DFT scans (HCN 11, NH3 9, H2CO 11 points; the CH4 C-H stretch has 9 rigid points), so the paths are the same in both notebooks;
* GFN2-xTB minimum, 21-point relaxed scan and peak (scout geometries) for every reaction;
* a per-structure SHA-256 over the canonical coordinates (rounded to 1e-8 angstrom), a per-path SHA-256, and one SHA-256 over the whole pack.

## Levels

| item | level |
|---|---|
| DFT geometries, Hessians, scans | DF-B3LYP-D3(BJ), integration grid level 3, SCF conv_tol 1e-9; def2-TZVP, def2-TZVPD for the SN2 species (`basis` in every structure) |
| scout geometries | GFN2-xTB (tblite), geomeTRIC relaxed scan along the same coordinate |
| gradient test on Colab (cell C0) | the same DF-B3LYP-D3(BJ) energy and gradient; pass when RMS gradient < 3e-4 hartree/bohr (the builder's own target is 1.5e-4) |

## Versions and command

Hardware: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 cores, 16 GB RAM  
Versions: `pyscf=2.15.0;geometric=1.1.1;tblite=0.8.0;ase=3.29.0;numpy=2.5.4;scipy=1.18.1;pandas=3.1.1;pyscf-dispersion=1.7.0;python=3.13.18`

```
python build_geompack.py --work <dir> --stage core --out geompack/b1v5_geompack_1h.json     # six reactions + CH4
python build_geompack.py --work <dir> --stage ext78 --out geompack/b1v5_geompack_1h.json    # adds F- + CH3Cl and OH- + CH3F (same work folder)
python build_geompack.py --work <dir> --assemble --out geompack/b1v5_geompack_1h.json       # write the pack file again, no calculation
```

The build is resumable (`<dir>/pack_work.json`). The notebooks embed the pack as gzip + base64 in one cell with the constant `GEOMPACK_SHA256`; cell 6 refuses a pack whose hash differs. `REBUILD_GEOMS = True` (CPU notebook) rebuilds everything on Colab with the same code (about 1.5 to 2 hours there, not part of the hour).

## Structures

| system | species | level | kind | basis | imaginary modes | lowest freq (cm-1) | RMS grad (hartree/bohr) | flags |
|---|---|---|---|---|---|---|---|---|
| ch3nc | CH3CN | own_b3lyp | min | def2-TZVP | 0 | 383.8 | 3.0e-05 |  |
| ch3nc | CH3NC | own_b3lyp | min | def2-TZVP | 0 | 276.4 | 3.9e-05 |  |
| ch3nc | TS | own_b3lyp | ts | def2-TZVP | 1 | -442.3 | 4.3e-05 |  |
| ch4 | CH4 | own_b3lyp | min | def2-TZVP | 0 | 1341.3 | 1.0e-07 |  |
| h2co | CO | own_b3lyp | min | def2-TZVP | 0 | 2213.9 | 8.4e-06 |  |
| h2co | H2 | own_b3lyp | min | def2-TZVP | 0 | 4415.3 | 3.9e-07 |  |
| h2co | H2CO | own_b3lyp | min | def2-TZVP | 0 | 1202.4 | 6.7e-06 |  |
| h2co | TS | own_b3lyp | ts | def2-TZVP | 1 | -1869.3 | 6.3e-06 |  |
| hcn | HCN | own_b3lyp | min | def2-TZVP | 0 | 766.8 | 1.4e-05 |  |
| hcn | HNC | own_b3lyp | min | def2-TZVP | 0 | 461.3 | 6.0e-06 |  |
| hcn | TS | own_b3lyp | ts | def2-TZVP | 1 | -1115.4 | 8.8e-06 |  |
| nh3 | NH3 | own_b3lyp | min | def2-TZVP | 0 | 1028.1 | 8.1e-05 |  |
| nh3 | TS | own_b3lyp | ts | def2-TZVP | 1 | -798.7 | 8.5e-07 |  |
| sn2_clcl | CH3Cl | own_b3lyp | min | def2-TZVPD | 0 | 699.1 | 5.8e-06 |  |
| sn2_clcl | Cl- | own_b3lyp | min | def2-TZVPD | 0 | None | 4.2e-15 |  |
| sn2_clcl | RC | own_b3lyp | min | def2-TZVPD | 0 | 73.9 | 5.7e-06 | Hessian at grid level 5; extra-mode retry done |
| sn2_clcl | TS | own_b3lyp | ts | def2-TZVPD | 1 | -362.0 | 1.3e-04 | Hessian at grid level 5; extra-mode retry done |
| sn2_fcl | CH3Cl | own_b3lyp | min | def2-TZVPD | 0 | 699.1 | 5.8e-06 |  |
| sn2_fcl | CH3F | own_b3lyp | min | def2-TZVPD | 0 | 1045.2 | 1.9e-06 |  |
| sn2_fcl | Cl- | own_b3lyp | min | def2-TZVPD | 0 | None | 4.6e-15 |  |
| sn2_fcl | F- | own_b3lyp | min | def2-TZVPD | 0 | None | 7.5e-16 |  |
| sn2_fcl | PC | own_b3lyp | min | def2-TZVPD | 0 | 88.0 | 6.9e-05 | Hessian at grid level 5; extra-mode retry done |
| sn2_fcl | RC | own_b3lyp | min | def2-TZVPD | 0 | 97.8 | 1.2e-05 |  |
| sn2_fcl | TS | own_b3lyp | ts | def2-TZVPD | 1 | -132.7 | 7.8e-06 |  |
| sn2_ff | CH3F | own_b3lyp | min | def2-TZVPD | 0 | 1045.2 | 1.9e-06 |  |
| sn2_ff | F- | own_b3lyp | min | def2-TZVPD | 0 | None | 2.8e-16 |  |
| sn2_ff | RC | own_b3lyp | min | def2-TZVPD | 0 | 116.0 | 5.2e-05 |  |
| sn2_ff | TS | own_b3lyp | ts | def2-TZVPD | 1 | -475.7 | 4.9e-07 |  |
| sn2_ohf | CH3F | own_b3lyp | min | def2-TZVPD | 0 | 1045.2 | 1.9e-06 |  |
| sn2_ohf | CH3OH | own_b3lyp | min | def2-TZVPD | 0 | 291.6 | 1.4e-05 | extra-mode retry done |
| sn2_ohf | F- | own_b3lyp | min | def2-TZVPD | 0 | None | 2.2e-16 |  |
| sn2_ohf | OH- | own_b3lyp | min | def2-TZVPD | 0 | 3738.0 | 6.3e-06 |  |
| sn2_ohf | PC | own_b3lyp | min | def2-TZVPD | 0 | 72.7 | 2.2e-05 | extra-mode retry done |
| sn2_ohf | RC | own_b3lyp | min | def2-TZVPD | 0 | 76.1 | 2.3e-05 | extra-mode retry done |
| sn2_ohf | TS | own_b3lyp | ts | def2-TZVPD | 1 | -498.1 | 1.5e-05 | extra-mode retry done |

## Paths

| system | label | points | level | basis |
|---|---|---|---|---|
| ch4 | C-H stretch (A) | 9 | own_scan | def2-TZVP |
| h2co | H...H distance (A) | 11 | own_scan | def2-TZVP |
| hcn | H-C-N angle (deg) | 11 | own_scan | def2-TZVP |
| nh3 | N-H1-H2-H3 dihedral (deg) | 9 | own_scan | def2-TZVP |

## xTB scout sets

| system | scan points | peak index | peak found |
|---|---|---|---|
| ch3nc | 2 | None | yes |
| h2co | 2 | None | no |
| hcn | 2 | None | yes |
| nh3 | 2 | None | yes |
| sn2_clcl | 2 | None | yes |
| sn2_fcl | 2 | None | no |
| sn2_ff | 2 | None | yes |
| sn2_ohf | 2 | None | yes |

## Measured times of every build step (sandbox, 4 cores, seconds)

Measured on the sandbox machine above; some steps ran while short test jobs shared the four cores (the numbers are upper bounds for a quiet machine). Colab (2 vCPU) is expected to take about twice as long.

| system | species | step | seconds | note |
|---|---|---|---|---|
| hcn | HCN | optimisation | 6.6 | kind=min converged=True |
| hcn | HCN | frequencies + zero-point energy | 17.4 | imag=0 |
| hcn | HCN | gradient at the stored geometry | 2.8 | rms=1.44e-05 max=3.68e-05 |
| hcn | HNC | optimisation | 7.5 | kind=min converged=True |
| hcn | HNC | frequencies + zero-point energy | 10.5 | imag=0 |
| hcn | HNC | gradient at the stored geometry | 2.9 | rms=6.02e-06 max=1.30e-05 |
| hcn | TS | optimisation | 39.1 | kind=ts converged=True |
| hcn | TS | frequencies + zero-point energy | 35.6 | imag=1 |
| hcn | TS | gradient at the stored geometry | 2.9 | rms=8.81e-06 max=1.80e-05 |
| hcn | path | relaxed DFT scan | 360.4 | 11 points, peak index 6 |
| nh3 | NH3 | optimisation | 5.0 | kind=min converged=True |
| nh3 | NH3 | frequencies + zero-point energy | 12.4 | imag=0 |
| nh3 | NH3 | gradient at the stored geometry | 1.9 | rms=8.06e-05 max=2.04e-04 |
| nh3 | TS | optimisation | 11.8 | kind=ts converged=True |
| nh3 | TS | frequencies + zero-point energy | 8.6 | imag=1 |
| nh3 | TS | gradient at the stored geometry | 1.8 | rms=8.52e-07 max=1.99e-06 |
| nh3 | path | relaxed DFT scan | 159.6 | 9 points, peak index 4 |
| h2co | H2CO | optimisation | 15.8 | kind=min converged=True |
| h2co | H2CO | frequencies + zero-point energy | 29.6 | imag=0 |
| h2co | H2CO | gradient at the stored geometry | 4.0 | rms=6.70e-06 max=1.71e-05 |
| h2co | TS | optimisation | 74.8 | kind=ts converged=True |
| h2co | TS | frequencies + zero-point energy | 58.1 | imag=1 |
| h2co | TS | gradient at the stored geometry | 5.6 | rms=6.33e-06 max=1.68e-05 |
| h2co | H2 | optimisation | 6.1 | kind=min converged=True |
| h2co | H2 | frequencies + zero-point energy | 6.8 | imag=0 |
| h2co | H2 | gradient at the stored geometry | 2.4 | rms=3.94e-07 max=6.83e-07 |
| h2co | CO | optimisation | 10.8 | kind=min converged=True |
| h2co | CO | frequencies + zero-point energy | 25.7 | imag=0 |
| h2co | CO | gradient at the stored geometry | 4.1 | rms=8.35e-06 max=1.86e-05 |
| h2co | path | relaxed DFT scan | 595.0 | 11 points, peak index 5 |
| ch3nc | CH3NC | optimisation | 26.1 | kind=min converged=True |
| ch3nc | CH3NC | frequencies + zero-point energy | 100.9 | imag=0 |
| ch3nc | CH3NC | gradient at the stored geometry | 7.8 | rms=3.87e-05 max=1.21e-04 |
| ch3nc | CH3CN | optimisation | 27.2 | kind=min converged=True |
| ch3nc | CH3CN | frequencies + zero-point energy | 81.6 | imag=0 |
| ch3nc | CH3CN | gradient at the stored geometry | 7.9 | rms=2.99e-05 max=9.36e-05 |
| ch3nc | TS | optimisation | 296.8 | kind=ts converged=True |
| ch3nc | TS | frequencies + zero-point energy | 69.2 | imag=1 |
| ch3nc | TS | gradient at the stored geometry | 8.3 | rms=4.29e-05 max=1.03e-04 |
| sn2_ff | F- | optimisation | 0.0 | kind=min converged=True |
| sn2_ff | F- | frequencies + zero-point energy | 1.4 | imag=0 |
| sn2_ff | F- | gradient at the stored geometry | 1.7 | rms=2.82e-16 max=4.10e-16 |
| sn2_ff | CH3F | optimisation | 24.7 | kind=min converged=True |
| sn2_ff | CH3F | frequencies + zero-point energy | 45.7 | imag=0 |
| sn2_ff | CH3F | gradient at the stored geometry | 6.6 | rms=1.94e-06 max=5.75e-06 |
| sn2_ff | RC | optimisation | 51.8 | kind=min converged=True |
| sn2_ff | RC | frequencies + zero-point energy | 108.0 | imag=0 |
| sn2_ff | RC | gradient at the stored geometry | 9.3 | rms=5.24e-05 max=1.22e-04 |
| sn2_ff | TS | optimisation | 140.6 | kind=ts converged=True |
| sn2_ff | TS | frequencies + zero-point energy | 127.5 | imag=1 |
| sn2_ff | TS | gradient at the stored geometry | 8.4 | rms=4.89e-07 max=1.36e-06 |
| sn2_clcl | Cl- | optimisation | 0.0 | kind=min converged=True |
| sn2_clcl | Cl- | frequencies + zero-point energy | 1.4 | imag=0 |
| sn2_clcl | Cl- | gradient at the stored geometry | 1.7 | rms=4.24e-15 max=5.54e-15 |
| sn2_clcl | CH3Cl | optimisation | 21.9 | kind=min converged=True |
| sn2_clcl | CH3Cl | frequencies + zero-point energy | 104.8 | imag=0 |
| sn2_clcl | CH3Cl | gradient at the stored geometry | 7.0 | rms=5.83e-06 max=1.27e-05 |
| sn2_clcl | RC | optimisation | 58.1 | kind=min converged=True |
| sn2_clcl | RC | frequencies + zero-point energy (+ retry along the extra mode) | 574.9 | imag=3 |
| sn2_clcl | RC | gradient at the stored geometry | 12.9 | rms=5.69e-06 max=1.52e-05 |
| sn2_clcl | TS | optimisation | 453.0 | kind=ts converged=True |
| sn2_clcl | TS | frequencies + zero-point energy (+ retry along the extra mode) | 1234.8 | imag=2 |
| sn2_clcl | TS | gradient at the stored geometry | 12.0 | rms=1.33e-04 max=3.99e-04 |
| sn2_clcl | TS | scan-peak fallback for the TS | 2265.3 | point 6 |
| ch4 | CH4 | optimisation | 4.8 | kind=min converged=True |
| ch4 | CH4 | frequencies + zero-point energy | 10.1 | imag=0 |
| ch4 | CH4 | gradient at the stored geometry | 2.1 | rms=1.01e-07 max=1.13e-07 |
| hcn | HCN | xTB minimum | 0.1 | converged=True |
| hcn | TS | xTB relaxed scan + peak | 1.7 | 21 points, peak 13 |
| nh3 | NH3 | xTB minimum | 0.0 | converged=True |
| nh3 | TS | xTB relaxed scan + peak | 1.3 | 21 points, peak 10 |
| h2co | H2CO | xTB minimum | 0.0 | converged=True |
| h2co | TS | xTB relaxed scan + peak | 1.3 | 21 points, peak None |
| ch3nc | CH3NC | xTB minimum | 0.0 | converged=True |
| ch3nc | TS | xTB relaxed scan + peak | 3.2 | 21 points, peak 14 |
| sn2_ff | RC | xTB minimum | 0.1 | converged=True |
| sn2_ff | TS | xTB relaxed scan + peak | 1.7 | 21 points, peak 14 |
| sn2_clcl | RC | xTB minimum | 0.1 | converged=True |
| sn2_clcl | TS | xTB relaxed scan + peak | 2.1 | 21 points, peak 13 |
| sn2_clcl | TS | Hessian repeated at grid level 5 | 744.1 | imag 2 -> 1 |
| sn2_clcl | RC | Hessian repeated at grid level 5 | 772.8 | imag 3 -> 0 |
| sn2_fcl | F- | optimisation | 0.0 | kind=min converged=True |
| sn2_fcl | F- | frequencies + zero-point energy | 1.4 | imag=0 |
| sn2_fcl | F- | gradient at the stored geometry | 1.6 | rms=7.48e-16 max=9.39e-16 |
| sn2_fcl | CH3Cl | optimisation | 21.4 | kind=min converged=True |
| sn2_fcl | CH3Cl | frequencies + zero-point energy | 92.7 | imag=0 |
| sn2_fcl | CH3Cl | gradient at the stored geometry | 6.8 | rms=5.83e-06 max=1.27e-05 |
| sn2_fcl | RC | optimisation | 88.4 | kind=min converged=True |
| sn2_fcl | RC | frequencies + zero-point energy | 218.2 | imag=0 |
| sn2_fcl | RC | gradient at the stored geometry | 10.7 | rms=1.20e-05 max=2.95e-05 |
| sn2_fcl | TS | optimisation | 428.3 | kind=ts converged=True |
| sn2_fcl | TS | frequencies + zero-point energy | 162.1 | imag=1 |
| sn2_fcl | TS | gradient at the stored geometry | 9.5 | rms=7.82e-06 max=2.04e-05 |
| sn2_fcl | PC | optimisation | 37.5 | kind=min converged=True |
| sn2_fcl | PC | frequencies + zero-point energy (+ retry along the extra mode) | 307.2 | imag=1 |
| sn2_fcl | PC | Hessian repeated at grid level 5 | 659.2 | imag 1 -> 0 |
| sn2_fcl | PC | gradient at the stored geometry | 9.1 | rms=6.86e-05 max=2.01e-04 |
| sn2_fcl | CH3F | optimisation | 23.6 | kind=min converged=True |
| sn2_fcl | CH3F | frequencies + zero-point energy | 34.6 | imag=0 |
| sn2_fcl | CH3F | gradient at the stored geometry | 5.9 | rms=1.94e-06 max=5.75e-06 |
| sn2_fcl | Cl- | optimisation | 0.0 | kind=min converged=True |
| sn2_fcl | Cl- | frequencies + zero-point energy | 1.3 | imag=0 |
| sn2_fcl | Cl- | gradient at the stored geometry | 1.6 | rms=4.56e-15 max=7.58e-15 |
| sn2_ohf | OH- | optimisation | 4.6 | kind=min converged=True |
| sn2_ohf | OH- | frequencies + zero-point energy | 3.8 | imag=0 |
| sn2_ohf | OH- | gradient at the stored geometry | 2.0 | rms=6.28e-06 max=1.16e-05 |
| sn2_ohf | CH3F | optimisation | 23.4 | kind=min converged=True |
| sn2_ohf | CH3F | frequencies + zero-point energy | 35.9 | imag=0 |
| sn2_ohf | CH3F | gradient at the stored geometry | 6.0 | rms=1.94e-06 max=5.75e-06 |
| sn2_ohf | RC | optimisation | 86.9 | kind=min converged=True |
| sn2_ohf | RC | frequencies + zero-point energy (+ retry along the extra mode) | 495.5 | imag=0 |
| sn2_ohf | RC | gradient at the stored geometry | 11.8 | rms=2.29e-05 max=4.22e-05 |
| sn2_ohf | TS | optimisation | 252.3 | kind=ts converged=True |
| sn2_ohf | TS | frequencies + zero-point energy (+ retry along the extra mode) | 526.5 | imag=1 |
| sn2_ohf | TS | gradient at the stored geometry | 10.9 | rms=1.52e-05 max=3.40e-05 |
| sn2_ohf | PC | optimisation | 73.4 | kind=min converged=True |
| sn2_ohf | PC | frequencies + zero-point energy (+ retry along the extra mode) | 545.8 | imag=0 |
| sn2_ohf | PC | gradient at the stored geometry | 13.3 | rms=2.17e-05 max=6.86e-05 |
| sn2_ohf | CH3OH | optimisation | 25.6 | kind=min converged=True |
| sn2_ohf | CH3OH | frequencies + zero-point energy (+ retry along the extra mode) | 317.8 | imag=0 |
| sn2_ohf | CH3OH | gradient at the stored geometry | 8.3 | rms=1.45e-05 max=4.53e-05 |
| sn2_ohf | F- | optimisation | 0.0 | kind=min converged=True |
| sn2_ohf | F- | frequencies + zero-point energy | 1.4 | imag=0 |
| sn2_ohf | F- | gradient at the stored geometry | 1.6 | rms=2.18e-16 max=2.88e-16 |
| sn2_fcl | RC | xTB minimum | 0.6 | converged=True |
| sn2_fcl | TS | xTB relaxed scan + peak | 0.7 | 21 points, peak None |
| sn2_ohf | RC | xTB minimum | 0.3 | converged=True |
| sn2_ohf | TS | xTB relaxed scan + peak | 25.8 | 21 points, peak 5 |

Sum of the recorded steps: 226 minutes (3.77 hours) in the sandbox, so about 7.5 hours on Colab for `REBUILD_GEOMS` (scan points are included in the steps above only where the builder times them as steps; see `pack_work.json` of the build for every scan point).

## Notes on this build

* **Cl- + CH3Cl, frequencies at integration grid level 5.** At the standard grid level 3 the Hessian of the TS had two imaginary modes (-373.5 and -134.1 cm-1) and the reactant complex three (-204.2, -68.2, -67.6 cm-1) even after the one retry along the extra modes. Repeating the Hessian at grid level 5 at the same stored geometries gave exactly one imaginary mode for the TS (-362.0 cm-1; the former -134 became +206 cm-1) and none for the complex. The extra modes were grid noise of very flat modes of an anion with diffuse functions. These two structures carry `hessian_grid: 5` (frequencies and zero-point energy from level 5; geometry, energy and gradient test stay at level 3). Structures without the key were analysed at level 3. The builder now does this repeat automatically (`PackBuilder.regrid_analysis`); the two structures here were patched with `build_geompack.py --regrid sn2_clcl:TS,sn2_clcl:RC` (744 s and 773 s in the sandbox).
* A scan-peak fallback for the Cl- + CH3Cl TS (11-point relaxed scan along the C-Cl distance, 2265 s, plus its Hessian) had been computed before the cause was found. It was dropped from the pack because the own TS has one imaginary mode at level 5; its times remain in the table below (they are not needed to rebuild the pack: a clean rebuild with the current code takes about 112 minutes in the sandbox instead of the 150 minutes recorded here).
* The TS gradient of Cl- + CH3Cl (RMS 1.3e-4 hartree/bohr) is the loosest of the pack; the builder's tightening threshold is 1.5e-4 and the Colab gradient test tolerance is 3e-4.
* The H2CO xTB scan has no peak (`no_peak`); its scout rows are skipped and the check (I12) says so.
* Timing caveat: the three heaviest build steps and the first sn2_clcl steps shared the four cores with short test jobs.
* **Reactions 7 and 8 (F- + CH3Cl, OH- + CH3F) are in the pack** (stage `ext78`, 76.6 minutes in the sandbox while the machine was otherwise idle apart from this build). Every structure has the wanted mode count (TS: one imaginary mode, minima: none). Three Hessians were repeated: the F- + CH3Cl ion-dipole product complex had one small imaginary mode at grid level 3 and none at grid level 5 (`hessian_grid: 5`); the OH- + CH3F TS, RC, PC and CH3OH needed the one retry along the extra mode (`retried: true`). Their gradient tests pass (worst RMS 6.9e-5 hartree/bohr).

