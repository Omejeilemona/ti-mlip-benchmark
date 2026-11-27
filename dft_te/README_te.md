# Finite-electron-temperature DFT campaign (Ti)

55 configurations, each computed at k_B*Te = 0.1, 0.5, 1.0, 1.5 and 2.0 eV: 275 GPAW jobs.
The 0.1 eV jobs are the reference for the electronic correction dF(R;Te) = F(R;Te) - F(R;0.1 eV),
so every Te uses the same plane waves, k-points and cells.

| Group | Configurations | Atoms |
|---|---|---|
| A  equation of state: hcp, bcc, fcc, omega at 0.94, 1.00, 1.06 V0 | 12 | 1-3 |
| B  hcp -> bcc Burgers path (strain + shuffle, 10 points) | 10 | 4 |
| C  hcp -> omega path | 8 | 1-3 |
| D  rattled omega | 4 | 6 |
| E  hot solid | 6 | 16 |
| F  liquid | 4 | 24 |
| G  surfaces | 3 | 16-20 |
| H  clusters in vacuum (Ti13 icosahedron, rattled, expanded; Ti19) | 4 | 13-19 |
| I  expanded hot solid (0.70 and 0.85 of normal density) | 4 | 16 |

Settings: PBE, 450 eV, k-point density 3.5 /A (1 k-point for clusters), Fermi-Dirac width = k_B*Te,
200% bands for kTe <= 1 eV and 300% for kTe >= 1.5 eV (checked: highest band empty). Each result stores
the Mermin free energy F (the energy the forces derive from), E0, forces, stress and the highest-band occupation.

## Run (from the repository folder, next to dft_v3)
    cd dft_te
    bash run_te.sh 3          # 3 workers; large jobs run one at a time
    bash status_te.sh         # progress
After a power cut or restart:
    cd dft_te
    python3 scripts/dft_worker_te.py --requeue
    bash run_te.sh 3
To stop all workers:  pkill -f dft_worker_te.py
The queue is already filled; scripts/gen_te_structures.py documents how it was built (needs ../dft_v3/done).
