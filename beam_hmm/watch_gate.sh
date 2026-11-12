#!/bin/bash
# Watchdog for a resumable gate run: relaunches it if it died (e.g. container restart) and exits
# with the report once all 40 per-fly results exist. Usage: watch_gate.sh <outdir>
cd "$(dirname "$0")"; OUT="$1"
while true; do
  n=$(ls "$OUT"/G[12]_fly*.json 2>/dev/null | wc -l)
  if [ "$n" -ge 40 ]; then python3 -I gate_d1.py "$OUT" --report; exit 0; fi
  if ! pgrep -f "gate_d1.py $OUT 4" > /dev/null; then
    echo "$(date -u) gate not running with $n/40 done: relaunching" >> "$OUT/watchdog.log"
    nohup python3 -I gate_d1.py "$OUT" 4 >> "$OUT/run.log" 2>&1 &
  fi
  sleep 120
done
