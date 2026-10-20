#!/bin/bash
# When each run finishes: append its raw log to the ledger (verbatim, marked AUTO), commit locally, push.
S=/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit; cd $S
done9=0; doneC=0
while [ $done9 = 0 ] || [ $doneC = 0 ]; do
  if [ $done9 = 0 ] && [ -f within_fly9.status ]; then
    { echo; echo "## $(date -u +%F) AUTO-APPENDED raw output — PRESPEC 9 (acc2ffd), within_fly9.py (interpretation to follow)"; sed 's/^/    /' within_fly9.log; } >> RESULTS_LEDGER.md
    git add -A; git -c user.name=audit -c user.email=audit@local commit -qm "Ledger: PRESPEC 9 raw (auto)"; ./sync_push.sh "PRESPEC 9 results (auto)"; done9=1; fi
  if [ $doneC = 0 ] && [ -f counts.status ]; then
    { echo; echo "## $(date -u +%F) AUTO-APPENDED raw output — PRESPEC 10 (ab86161), counts_analysis.py (interpretation to follow)"; sed 's/^/    /' counts.log; } >> RESULTS_LEDGER.md
    git add -A; git -c user.name=audit -c user.email=audit@local commit -qm "Ledger: PRESPEC 10 raw (auto)"; ./sync_push.sh "PRESPEC 10 results (auto)"; doneC=1; fi
  # interim push of partial logs every 30 min
  sleep 1800; [ $done9 = 0 ] || [ $doneC = 0 ] && ./sync_push.sh "interim logs (auto)" >/dev/null 2>&1
done
