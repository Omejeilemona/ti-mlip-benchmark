#!/bin/bash
# Constant-Te MD, cold lattice (300 K): kTe = 0 (control), 0.5, 1.0, 1.5, 2.0 eV; flexible (npt) and fixed (nvt) cell.
export OMP_NUM_THREADS=1
for m in npt nvt; do for t in 0.0 0.5 1.0 1.5 2.0; do echo "--kte $t --mode $m"; done; done | \
  xargs -P ${1:-3} -I{} bash -c 'python3 scripts/md_constant_te.py {} >> md_te_summary.txt 2>&1'
grep "kTe =" md_te_summary.txt | sort
