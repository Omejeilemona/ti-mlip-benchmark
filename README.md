# ti-mlip-benchmark

**Machine-learned interatomic potentials for titanium, including an electron-temperature-dependent correction for laser-excited Ti**

Ilemona S. Omeje, Moses O. Alhassan and Tatiana E. Itina
Surface Optimization Research & Consulting (SOREC) LLC, Hattiesburg, Mississippi, USA, and Laboratoire Hubert Curien, CNRS / Université Jean Monnet, Saint-Étienne, France

This repository contains the potentials, scripts and processed results of the paper

> I. S. Omeje, M. O. Alhassan and T. E. Itina, *Electron-temperature-dependent machine-learned potentials reveal hot-electron effects on phase transformations and nanoparticle fragmentation in titanium*, submitted to Physical Review Materials (2026).

The DFT data and an archived copy of this repository are available on Zenodo: **https://doi.org/10.5281/zenodo.24591569**

## Overview

The work has two parts.

1. **Ground-state potentials.** SNAP, POD and MACE potentials for Ti are fitted to a public PBE dataset and tested on 300 held-out configurations, on lattice, elastic, phase-stability, defect and surface properties, and on the melting point. The missing ω phase of the reference data is repaired with 93 targeted GPAW calculations (POD v3, MACE v3).
2. **Electron-temperature dependence.** 316 finite-electron-temperature DFT calculations (k_B T_e = 0.1–2 eV) give a correction ΔF(R; T_e) = F(R; T_e) − E₀(R), fitted as a linear SNAP model in 21 bins of 0.1 eV. The correction runs in standard LAMMPS through `pair_style hybrid/overlay` and is added to a ground-state potential. It is used in constant-T_e molecular dynamics of hcp Ti and in two-temperature-model MD (TTM-MD) of a laser-heated 2 nm Ti nanoparticle.

## Main results

| Quantity | Value |
|---|---|
| ω − hcp energy (DFT: −4 meV/atom) | POD v3: −3.7 meV/atom, MACE v3: −7.6 meV/atom |
| Melting point (experiment: 1941 K) | POD v3: 1945 ± 35 K, MACE v3: 1770 ± 20 K |
| bcc more stable than hcp (cold lattice) | above k_B T_e ≈ 0.3 eV |
| fcc becomes the lowest-energy phase | above k_B T_e ≈ 1.2 eV |
| Electronic pressure at 2 eV | about 40 GPa |
| hcp at 0.5 eV, lattice at 300 K, flexible cell | transforms within about 1 ps (fcc stacking, then bcc) |
| 2 nm nanoparticle, energy for complete disintegration | EAM: about 5 eV/atom, EAM + ΔF: about 3 eV/atom |

## Contents

| Folder | Contents |
|---|---|
| `potentials/snap` | SNAP (2J = 10) with ZBL core |
| `potentials/pod_v1`, `pod_v2_omega`, `pod_v3` | POD fitted to the reference data, + ω data, + all 93 new DFT configurations |
| `potentials/pod_v4` | POD v3 refitted with clusters and expanded cells; used with a ZBL core and a short-range wall (`wall.table`) |
| `potentials/mace`, `potentials/mace_v3` | Compact MACE (reference data) and MACE v3 (v3 data; also a TorchScript-compiled copy that loads on a CPU) |
| `potentials/te_snap` | Electron-temperature correction ΔF = F − E₀ (bins `Te00`–`Te20`); add to any ground-state potential |
| `potentials/te_snap_res` | Residual correction ΔF_res = F − E_POD v4, for use with POD v4 only |
| `dft_v3` | Ground-state GPAW campaign: structure generation, job queue and the 93 results |
| `dft_te` | Finite-T_e GPAW campaign: structure generation, workers and the results (`done/`) |
| `data_te` | Finite-T_e data collected per electron temperature, with `summary.csv` |
| `data_omega`, `data_mace` | ω data and the MACE training / validation / test sets |
| `scripts` | Data conversion, fitting, held-out tests, properties, melting, phase stability, stacking faults, constant-T_e MD |
| `nanoparticle` | 0D TTM-MD of a free Ti nanoparticle (`ttm0d.py`), Lin–Zhigilei C_e and G tables, results and figure |
| `md_te` | Constant-T_e MD results |
| `results/mace_v3` | MACE v3 held-out errors, properties, melting runs and training curves |
| `notebooks` | Google Colab notebooks for training and analysing MACE v3 |
| `inputs`, `results`, `figures` | Fit inputs, error tables, property files and figures |

## Using the potentials in LAMMPS

POD v3 (ground state):

```
pair_style pod
pair_coeff * * potentials/pod_v3/Ti_param.pod potentials/pod_v3/Ti_coefficients.pod Ti
```

POD v4 with the ZBL core and short-range wall (as in `scripts/fit_pod_v4.py`):

```
pair_style hybrid/overlay pod zbl 1.2 1.6 table linear 1401
pair_coeff * * pod potentials/pod_v4/Ti_param.pod potentials/pod_v4/Ti_coefficients.pod Ti
pair_coeff * * zbl 22 22
pair_coeff * * table potentials/pod_v4/wall.table WALL 1.6
```

EAM + electron-temperature correction at k_B T_e = 0.5 eV (bin `Te05`). The EAM file (Mendelev *et al.* 2016) is not included; download `Ti1_Mendelev2016.eam.fs` from the [NIST Interatomic Potentials Repository](https://www.ctcms.nist.gov/potentials/) into `potentials/`.

```
pair_style hybrid/overlay eam/fs snap
pair_coeff * * eam/fs potentials/Ti1_Mendelev2016.eam.fs Ti
pair_coeff * * snap potentials/te_snap/Ti_te.snapcoeff potentials/te_snap/Ti_te.snapparam Te05
```

The bin index is 10 × k_B T_e in eV (`Te00` to `Te20`). In TTM-MD the bin is updated as T_e changes; see `nanoparticle/ttm0d.py`.

MACE (Python, through ASE):

```python
from mace.calculators import MACECalculator
atoms.calc = MACECalculator(model_paths='potentials/mace_v3/ti_mace_v3_stagetwo.model', device='cpu')
```

## Reproducing the results

The reference DFT data (Wen *et al.* 2021, as curated by Röcken and Zavadlav 2024, Apache-2.0) is not redistributed here; `./get_data.sh` downloads it. All scripts run from the repository root unless stated otherwise.

| Step | Command |
|---|---|
| Reference data and held-out test set | `./get_data.sh`, `python scripts/convert.py ...`, `python scripts/make_testset.py` |
| SNAP and POD v1 | `python -m fitsnap3 inputs/Ti-snap.in`, `lmp -in inputs/in.fitpod` |
| v3 DFT campaign and POD v3 | `cd dft_v3 && bash run_dft.sh 2`, then `python scripts/fit_v3.py` |
| Held-out tests and properties | `python scripts/heldout.py ...`, `python scripts/props.py ...` |
| Melting point (POD v3) | `python scripts/melt_pod.py ...`, `python scripts/q6.py ...` |
| MACE v3 | `python scripts/mace_prepare_v3.py`, then `notebooks/Ti_MACE_v3_colab.ipynb` on a GPU |
| Finite-T_e DFT campaign | `cd dft_te && bash run_te.sh 3`, then `python scripts/collect_te.py` |
| Phase stability, density of states, C_e | `python scripts/te_phase_analysis.py`, `python scripts/dos_phases.py`, `python scripts/ce_dft.py` |
| POD v4 and the T_e corrections | `python scripts/fit_pod_v4.py`, `python scripts/fit_te_snap.py`, `python scripts/fit_te_snap_res.py` |
| Constant-T_e MD | `bash run_constant_te_res.sh 3` |
| Nanoparticle TTM-MD | `cd nanoparticle && bash run_all.sh 3 && python fig_np_compare.py 3.5` |

Tested with LAMMPS 22 Jul 2025 (ML-POD, ML-SNAP, MANYBODY, EXTRA-COMPUTE), FitSNAP, GPAW 26.7, ASE and mace-torch.

## Citation

If you use these potentials or data, please cite the paper above and the Zenodo record. GitHub's "Cite this repository" button uses `CITATION.cff`.

## License

The code is released under the MIT License (`LICENSE`). The new DFT data and the potentials archived on Zenodo are released under CC BY 4.0. The configurations in `data_mace/` are derived from the reference dataset and remain under its Apache-2.0 license (Röcken and Zavadlav 2024; Wen *et al.* 2021).
