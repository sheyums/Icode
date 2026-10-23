#!/bin/bash
# Remaining PRESPEC 20b checks, then PRESPEC 20c, then theta sensitivity. Each stage logs to its own .out.
cd "$(dirname "$0")"
python3 sqw2.py boot > boot2.out 2>&1
python3 sqw2.py split > split2.out 2>&1
python3 sqw2.py adequacy > adeq2.out 2>&1
echo checks_done > stage.status
python3 sqw3.py synth > sqw3_synth.out 2>&1
python3 sqw3.py realistic > sqw3_realistic.out 2>&1
python3 sqw3.py fit 0.5 > fit3.out 2>&1
python3 sqw3.py route 0.5 > route3.out 2>&1
echo 20c_done > stage.status
for th in 0.33 0.75; do python3 sqw2.py fit $th K3only > fit2_$th.out 2>&1; python3 sqw2.py route $th >> fit2_$th.out 2>&1; done
echo all_done > stage.status
