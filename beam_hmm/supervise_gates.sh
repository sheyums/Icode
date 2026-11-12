#!/bin/bash
# Keeps the queued gates running in order: D1-B (40 flies) then T1 (60 flies). Relaunches a gate that died.
# Runs under nohup so it outlives tool time limits; after a container restart the scheduled check-ins restart it.
cd "$(dirname "$0")"
run_until() {  # <script> <outdir> <n_expected>
  while [ "$(ls "$2"/*_fly*.json 2>/dev/null | wc -l)" -lt "$3" ]; do
    if ! pgrep -f "$1 $2 4" > /dev/null; then
      echo "$(date -u) launching $1 $2 ($(ls "$2"/*_fly*.json 2>/dev/null | wc -l)/$3 done)" >> supervisor.log
      nohup python3 -I "$1" "$2" 4 >> "$2/run.log" 2>&1 &
    fi
    sleep 120
  done
  python3 -I "$1" "$2" --report > "$2/REPORT.txt" 2>&1
  echo "$(date -u) $1 complete" >> supervisor.log
}
mkdir -p results/t1_gate
run_until gate_d1.py results/d1b_gate 40
run_until gate_t1.py results/t1_gate 60
echo "$(date -u) all gates complete" >> supervisor.log
