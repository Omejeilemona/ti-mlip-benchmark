#!/bin/bash
# Prints mean temperature and solid fraction per slab for each run.
cd melt_v2
for f in run_*_final.dump; do
  p=${f%_final.dump}
  grep "mean T" $p.out || true
  python3 ../scripts/q6.py $f | tail -n 1
  echo
done
echo "Valid runs: some slabs mostly solid (>0.6) AND some mostly liquid (<0.3)."
echo "The melting point is the mean T of the valid runs."
