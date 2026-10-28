#!/bin/bash
# usage: run_variant.sh name twojmax rcut quad eweight
n=$1; VW=${6:-1.0E-6}; mkdir -p v_$n; sed -e "s/^twojmax = .*/twojmax = $2/" -e "s/^rcutfac = .*/rcutfac = $3/" -e "s/^quadraticflag = .*/quadraticflag = $4/" \
 -e "s/^Init = .*/Init = 0.8 0.2 $5 1 ${VW}/" -e "s/^Bulk = .*/Bulk = 0.8 0.2 $5 1 0.0/" \
 -e "s#^dataPath = .*#dataPath = ../XYZ#" -e "s/^metrics = .*/metrics = metrics.md/" -e "s/^potential = .*/potential = Ti_snap/" Ti-snap.in > v_$n/in.in
cd v_$n && PYTHONPATH=${FITSNAP_DIR:-../FitSNAP} python3 -m fitsnap3 in.in > log 2>&1
grep -E "'\*ALL', 'Unweighted', 'Testing'" metrics.md
