#!/bin/bash
# Update 1 of the finite-Te campaign. Run from dft_te/:  bash apply_update.sh
set -e
pkill -f dft_worker_te.py || true; sleep 3
python3 scripts/dft_worker_te.py --requeue
# 1. small clusters (group J) at kTe = 0.1, 0.5, 1.0 eV
mv queue_add/*.xyz queue/
# 2. clusters (H) and surfaces (G) are not computed at kTe = 1.5 and 2.0 eV: electrons are thermally emitted
#    into the vacuum, which a finite DFT box cannot represent
mkdir -p skipped
for f in queue/Te1.5__*__[GH]_*.xyz queue/Te2.0__*__[GH]_*.xyz; do [ -e "$f" ] && mv "$f" skipped/; done
echo "queue: $(ls queue | wc -l) jobs   skipped: $(ls skipped | wc -l)"
bash run_te.sh 3
