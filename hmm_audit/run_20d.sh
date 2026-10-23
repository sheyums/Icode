#!/bin/bash
cd "$(dirname "$0")"
until grep -q all_done stage.status 2>/dev/null; do sleep 60; done
python3 sqw4.py a > sqw4_a.out 2>&1
python3 sqw4.py b > sqw4_b.out 2>&1
echo 20d_done > stage_d.status
python3 sqw5.py > sqw5.out 2>&1
echo 20e_null_done > stage_e.status
