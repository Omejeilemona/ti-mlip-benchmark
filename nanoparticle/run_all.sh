#!/bin/bash
# Nanoparticle TTM-MD: EAM and EAM + DFT Te correction, 7 absorbed energies, 2 nm particle, 30 ps. Run from nanoparticle/.
N=${1:-3}; export OMP_NUM_THREADS=1; jobs=()
for pot in eam eamte; do for e in 1.5 2 2.5 3 3.5 4 5; do
  jobs+=("python3 ttm0d.py --pot $pot --eabs $e --tend 30 --dumpat 0,1,2,4,8,15,30 --out ${pot}_d2_e${e}")
done; done
printf '%s\n' "${jobs[@]}" | xargs -P $N -I{} bash -c '{} >> run_all.log 2>&1'
grep -h "atoms," run_all.log | sort
