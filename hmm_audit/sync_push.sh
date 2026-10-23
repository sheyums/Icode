#!/bin/bash
# Copy audit files (no data, no .npy) into the repo and push to the session branch.
S=/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit
cd /home/user/Icode || exit 1
for f in $S/*; do case "$f" in *.npy|*/__pycache__|*.status) ;; *) [ -f "$f" ] && cp "$f" hmm_audit/ ;; esac; done
for d in golden figures; do [ -d $S/$d ] && mkdir -p hmm_audit/$d && cp $S/$d/* hmm_audit/$d/; done
(cd $S && git log --format='%h %ci %s' --reverse) > hmm_audit/AUDIT_GIT_LOG.txt
git add hmm_audit
git diff --cached --quiet && exit 0
git commit -q -m "hmm_audit: sync ($1)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01RdsgANv31WkZVEr9KjfHyX"
for d in 2 4 8 16; do git push -u origin claude/trusting-allen-7nxrsf && exit 0; sleep $d; done
