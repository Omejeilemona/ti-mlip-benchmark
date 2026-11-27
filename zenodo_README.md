# Data and potentials: electron-temperature-dependent machine-learned potentials for titanium

Ilemona S. Omeje, Moses O. Alhassan and Tatiana E. Itina

This record contains the DFT data, interatomic potentials and processed simulation results of the paper

> I. S. Omeje, M. O. Alhassan and T. E. Itina, *Electron-temperature-dependent machine-learned potentials reveal hot-electron effects on phase transformations and nanoparticle fragmentation in titanium*, submitted to Physical Review Materials (2026).

Code: https://github.com/Omejeilemona/ti-mlip-benchmark (a copy of release v1.0 is included as `code_v1.0.zip`).

## Contents

| Folder | Contents |
|---|---|
| `potentials/` | All potentials of the paper: SNAP; POD v1, v2 (ω), v3, v4 (+ `wall.table`); compact MACE and MACE v3; the electron-temperature correction `te_snap` (ΔF = F − E₀, bins Te00–Te20) and the residual correction `te_snap_res` (ΔF_res = F − E_POD v4) |
| `dft_groundstate_v3/` | 93 GPAW-PBE ground-state configurations (ω phase, hcp–ω path, point defects, liquid, hot solid, surfaces, alignment cells), extended XYZ with energy, forces and stress |
| `dft_omega/` | Earlier ω-phase and alignment calculations (POD v2) |
| `dft_finite_te/done/` | 414 finite-electron-temperature GPAW-PBE calculations, one extended XYZ file per configuration and electron temperature (see below) |
| `dft_finite_te/by_temperature/` | The main campaign collected per electron temperature (0.1, 0.5, 1.0, 1.5, 2.0 eV), with `summary.csv` |
| `dft_finite_te/logs/` and `dft_groundstate_v3/logs/` | GPAW output of every calculation |
| `mace_training_sets/` | Training, validation and test sets of the compact MACE (`compact/`) and MACE v3 (`v3/`) |
| `results/constant_te_md/` | Constant-T_e MD of hcp Ti (time series, CSV) |
| `results/nanoparticle/` | 0D TTM-MD of the 2 nm nanoparticle: time series (`ttm.csv`) and settings of every run, and the switching-energy decomposition |
| `results/melting_pod_v3/` | Temperature series of the solid–liquid coexistence runs |
| `results/properties/` | Phase energies vs T_e, densities of states, electron heat capacity, error tables and property files |
| `ttm_inputs/` | Electron heat capacity and electron–phonon coupling tables of Ti from Lin, Zhigilei and Celli (Phys. Rev. B 77, 075133, 2008) in the format read by `ttm0d.py` |

## Finite-electron-temperature data

File names follow `Te<kTe in eV>__<atoms>__<group>_<name>.xyz`, for example `Te0.5__016__E_hotsolid_2.xyz`. The header of each file gives the electron temperature (`kTe`, eV), the Mermin free energy `F` (also stored as `energy` and `free_energy`, so that the forces are its exact derivatives), the zero-width energy `E0`, the stress, the number of bands (`nbands`) and the occupation of the highest band (`top_band_occ`); the atom lines give positions and forces. Calculations used GPAW (PBE, plane waves 450 eV, Γ-centred k-point density 3.5 Å⁻¹, Fermi–Dirac occupations at the stated T_e, bands increased with T_e until the highest band is empty to better than 10⁻⁴).

| Set | k_B T_e (eV) | Calculations |
|---|---|---|
| Main campaign | 0.1, 0.5, 1.0, 1.5, 2.0 | 74 + 73 + 73 + 48 + 48 = 316 |
| hcp–bcc crossover (equation of state and Burgers path) | 0.2, 0.3, 0.4 | 3 × 22 = 66 |
| Basal generalised stacking fault | 0.1, 0.5, 1.0 | 11 + 11 + 10 = 32 |
| Total | | 414 |

Groups: A equation of state (hcp, bcc, fcc, ω), B hcp–bcc Burgers path, C hcp–ω path, D rattled ω, E hot solid, F liquid, G surfaces, H clusters of 13 and 19 atoms, I expanded hot solid (0.70 and 0.85 of normal density), J atom, dimer, trimers, tetramers and clusters of 5–8 atoms, K basal stacking fault.

## Units and conventions

Energies in eV, lengths in Å, forces in eV/Å, stress in eV/Å³ (ASE convention). Ground-state GPAW energies in `dft_groundstate_v3/` are raw GPAW values; to combine them with the reference VASP data, shift by Δ = −1.194 eV/atom (see `scripts/fit_v3.py`).

## Reference data not included

The 5,704 reference configurations (Wen *et al.*, npj Comput. Mater. 7, 206, 2021; curated by Röcken and Zavadlav, npj Comput. Mater. 10, 69, 2024) are available from https://github.com/tummfm/Fused-EXP-DFT-MLP under Apache-2.0. The MACE training sets in this record contain configurations from that dataset and remain under Apache-2.0. The EAM potential Ti1 of Mendelev *et al.* (2016) is available from the NIST Interatomic Potentials Repository.

## License

CC BY 4.0, except `mace_training_sets/` (Apache-2.0, see above) and `code_v1.0.zip` (MIT).
