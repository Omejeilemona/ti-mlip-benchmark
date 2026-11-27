#!/bin/bash
echo "queue: $(ls queue | wc -l)   running: $(ls running/*.xyz 2>/dev/null | wc -l)   done: $(ls done | wc -l)   failed: $(ls failed | wc -l)"
echo "done per Te:"; for t in 0.1 0.5 1.0 1.5 2.0; do
  tot=$(( $(ls done | grep -c "^Te${t}__") + $(ls queue | grep -c "^Te${t}__") + $(ls running | grep -c "__Te${t}__") ))
  echo "  kTe=$t eV: $(ls done | grep -c "^Te${t}__")/$tot"; done
echo "workers alive: $(pgrep -f dft_worker_te.py | wc -l)"
echo "last results:"; tail -q -n 2 logs/worker_*.out 2>/dev/null | grep -v '^$' | tail -6
grep -h WARNING logs/worker_*.out 2>/dev/null | tail -3
