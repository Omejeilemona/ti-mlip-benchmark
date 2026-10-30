#!/bin/bash
# Progress of the DFT queue.
q=$(ls queue/*.xyz 2>/dev/null | wc -l); r=$(ls running/* 2>/dev/null | wc -l)
d=$(ls done/*.xyz 2>/dev/null | wc -l); f=$(ls failed/* 2>/dev/null | wc -l)
echo "waiting $q | running $r | done $d | failed $f | workers alive: $(pgrep -fc dft_worker.py)"
echo "--- done per category:"
ls done/*.xyz 2>/dev/null | xargs -n1 basename | cut -d_ -f2 | sort | uniq -c
echo "--- last finished:"
cat logs/worker_*.out 2>/dev/null | grep "E/atom" | sort | tail -n 5
