#!/bin/bash
cd ~/ti-mlip-benchmark/nanoparticle
for d in eam*_d2_e*; do printf "%-14s %6s ps\n" $d "$(tail -1 $d/ttm.csv 2>/dev/null | cut -d, -f1)"; done
echo; echo "running: $(pgrep -fc ttm0d.py)   finished: $(grep -c 'atoms,' run_all.log) / 14"
