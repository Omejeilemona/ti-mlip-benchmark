#!/bin/bash
# Start N DFT workers in the background (default 4). Each uses one CPU core.
N=${1:-4}
export OMP_NUM_THREADS=1
mkdir -p logs
for i in $(seq 1 $N); do
  nohup /usr/bin/python3 scripts/dft_worker.py > logs/worker_$i.out 2>&1 &
done
echo "Started $N workers. Check with: bash status.sh"
