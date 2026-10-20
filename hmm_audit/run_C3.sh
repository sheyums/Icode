#!/bin/bash
cd /tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit
python3 -I -W ignore counts_C3_resume.py > counts_C3.log 2>&1 && echo done > counts_C3.status
{ echo; echo "## $(date -u +%F) AUTO-APPENDED raw output — PRESPEC 10 C3 (counts_C3_resume.py), interpretation to follow"; sed 's/^/    /' counts_C3.log; } >> RESULTS_LEDGER.md
git add -A; git -c user.name=audit -c user.email=audit@local commit -qm "Ledger: PRESPEC 10 C3 raw (auto)"; ./sync_push.sh "PRESPEC 10 C3 results (auto)"
