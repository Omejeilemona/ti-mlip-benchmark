#!/bin/bash
# Constant-Te MD with the residual correction (POD v4 + te_snap_res), cold lattice 300 K
export OMP_NUM_THREADS=1
printf '%s\n' "--kte 0.5 --mode nvt" "--kte 1.0 --mode nvt" "--kte 0.5 --mode npt" | \
  xargs -P 3 -I{} bash -c 'python3 scripts/md_constant_te.py {} --corr potentials/te_snap_res >> md_te_res_summary.txt 2>&1'
grep "kTe =" md_te_res_summary.txt | sort
