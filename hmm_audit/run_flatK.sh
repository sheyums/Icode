#!/bin/bash
cd /tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/hmm_audit
python3 -I -W ignore flatK.py >> flatK.log 2>&1 && echo done > flatK.status
