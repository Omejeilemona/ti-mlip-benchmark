#!/bin/bash
# Start N workers in the background (default 3). Each uses one CPU core.
if pgrep -f dft_worker_te.py >/dev/null; then echo "Workers already running ($(pgrep -f dft_worker_te.py | wc -l)). Not starting more."; exit 1; fi
N=${1:-3}
export OMP_NUM_THREADS=1
mkdir -p logs
for i in $(seq 1 $N); do
  nohup python3 scripts/dft_worker_te.py > logs/worker_$i.out 2>&1 &
done
echo "Started $N workers. Check progress with: bash status_te.sh"
