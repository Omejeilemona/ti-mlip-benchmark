#!/bin/bash
# Four independent POD v2 coexistence runs in parallel (one CPU core each).
# A run is valid only if the final snapshot still contains BOTH solid and liquid (check with q6.py).
set -e
POT=potentials/pod_v2_omega
mkdir -p melt_v2 && cd melt_v2
for T in 1850 1900 1950 2000; do
  nohup python3 ../scripts/melt_long.py ../$POT run_$T $T > run_$T.out 2>&1 &
done
echo "Started 4 runs in melt_v2/. Check progress with:  tail -n 2 melt_v2/run_*.log"
echo "When all have finished (run_*.out contains a 'mean T' line), run:  bash summarise_melt.sh"
