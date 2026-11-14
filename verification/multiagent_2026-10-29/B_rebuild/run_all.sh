#!/bin/sh
# Sequential runs (n_jobs=3 inside each). Logs to out/log_*.txt
PY=$HOME/.venvs/ete-fam/Scripts/python.exe
cd "$(dirname "$0")"
mkdir -p out
for o in hba1c weight sbp ed_ip; do $PY -W ignore b_analysis.py $o > out/log_$o.txt 2>&1; done
for o in hba1c weight sbp ed_ip; do $PY -W ignore b_analysis.py $o --pc > out/log_${o}_pc.txt 2>&1; done
