#!/bin/bash
# Reproduce recorded PKM msolve degrees. Prints: name run maxdeg rounds rows cols
R=/home/user/aburan28/crypto
SYS=$R/research/ic_next10_20261014/03-pkm-degree/systems
OUT=/tmp/claude-0/-home-user-ecbench/237dd0ae-3c31-5963-b1b5-f5f7c04af389/scratchpad/harness/out/pkm
for name in pkm_kummer_p2013265921_m2_t6_x464833983 \
            pkm_kummer_p2013265921_m3_t3_x388393539 \
            pkm_kummer_p2013265921_m4_t2_x92234565 \
            pkm_kummer_p2013265921_m3_t4_x334226296; do
  for run in 1 2; do
    log=$OUT/${name}_run${run}.log
    msolve -t 1 -v 2 -f $SYS/${name}.ms -o $OUT/${name}_run${run}.out > $log 2>&1
    # parse: data rows begin with an integer (the deg); take max of col1, count rounds, max rows/cols
    awk '/^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+ x [0-9]+/ {
           rounds++; if($1>D)D=$1; if($4>MR)MR=$4; c=$6; if(c>MC)MC=c }
         END{ printf "%s run%s maxdeg=%s rounds=%s rows=%s cols=%s\n", N, RUN, D, rounds, MR, MC }' N=$name RUN=$run $log
  done
done
