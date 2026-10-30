Ti v3 DFT campaign (GPAW, PBE, plane waves 450 eV, k-point density 5/A^-1, Fermi-Dirac 0.1 eV)
93 configurations: alignment 4, omega 15, transition paths 14, defects 15 (vacancies + interstitials),
liquid 20, hot solid 16, surfaces 9.

Run (from ~/ti-mlip-benchmark/dft_v3, environment active):
  bash run_dft.sh 3          # start 3 workers (one core each, ~1.5-3 GB RAM each)
  bash status.sh             # progress
  python3 scripts/dft_worker.py --requeue   # after a crash or reboot, put interrupted jobs back
  bash run_dft.sh 3          # then restart

Stop all workers:  pkill -f dft_worker.py   (then --requeue before restarting)
Results: done/*.xyz (energy, forces, stress). Send the done/ folder, or fit with
  cd ~/ti-mlip-benchmark && python3 scripts/fit_v3.py
