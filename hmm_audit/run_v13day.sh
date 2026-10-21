#!/bin/bash
cd /tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit
while [ ! -f validate13.status ]; do sleep 60; done
python3 -I -W ignore validate13_day.py >> validate13_day.log 2>&1 && echo done > validate13_day.status
git add -A; git -c user.name=audit -c user.email=audit@local commit -qm "validate13_day finished (auto)"; ./sync_push.sh "day-model validation (auto)"
