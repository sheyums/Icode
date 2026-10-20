#!/bin/bash
cd /tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit
while [ ! -f counts_C3.status ]; do sleep 1800; git add -A >/dev/null 2>&1; git -c user.name=audit -c user.email=audit@local commit -qm "C3 partial checkpoint (auto)" >/dev/null 2>&1; ./sync_push.sh "C3 partial checkpoint (auto)" >/dev/null 2>&1; done
