# MACE v3 results

Model: `potentials/mace_v3/ti_mace_v3_stagetwo.model` (trained on a GPU, 64 channels, L_max = 1, correlation order 3, two interactions, r_max = 5 Å).
`ti_mace_v3_stagetwo_compiled.model` is the same model compiled with TorchScript; it also loads on a CPU-only machine (see `notebooks/Ti_MACE_v3_analysis.ipynb`).

| File | Contents |
|---|---|
| `test_errors_v3.txt` | Held-out errors on the 300 test configurations (Init: 2.3 meV/atom, 0.006 eV/Å; Bulk: 6.5 meV/atom, 0.159 eV/Å) |
| `mace_v3_props.json` | Lattice, elastic constants, phase energies and vacancy formation energy |
| `melt_T*.dat` | Solid–liquid coexistence runs (time in ps, energy in eV/atom, temperature in K, volume in Å³/atom) |
| `melting_v3.txt` | Sign of dE/dt at each temperature |
| `training_log.txt`, `training_stage_*.png` | Training log and learning curves |

Melting point. At 1700 K the solid grows (dE/dt = −4.7 meV/atom/ps) and at 1800 K it melts (+0.82 meV/atom/ps); at 1775 K the slope is close to zero (+0.07). This brackets Tm = 1770 ± 20 K as reported in the paper. The 1900 K file is a shorter, preliminary 20 ps run whose slope (−0.08 meV/atom/ps) is within the noise of zero; it is inconclusive and was not used.
