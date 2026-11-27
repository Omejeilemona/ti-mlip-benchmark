#!/bin/bash
cd ~/ti-mlip-benchmark
echo "run              time(ps)  hcp    bcc    fcc    P(GPa)  V(A3/atom)"
for f in md_te/*.log; do
  c=${f%.log}.csv; n=$(basename ${f%.log})
  if [ -f $c ]; then tail -1 $c | awk -F, -v n=$n '{printf "%-16s %6s   %5s  %5s  %5s  %6s  %6s\n", n, $1, $2, $3, $4, $7, $8}'
  else echo "$n  equilibrating, step $(grep -E '^ +[0-9]+ ' $f | tail -1 | awk '{print $1}') of 5000"; fi
done
echo; echo "finished: old $(grep -c "kTe =" md_te_summary.txt 2>/dev/null)/10, new (_res) $(grep -c "kTe =" md_te_res_summary.txt 2>/dev/null)/3   $(date +%H:%M)"
