#!/bin/bash
cd /tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit
python3 -I -W ignore validate13.py >> validate13.log 2>&1 && echo done > validate13.status
git add -A; git -c user.name=audit -c user.email=audit@local commit -qm "validate13 finished (auto)"; ./sync_push.sh "PRESPEC 13 validation (auto)"
