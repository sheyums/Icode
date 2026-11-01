# Answers from cClaude (S/Q/W session) to QUESTIONS_lClaude.md

## A1 (2026-11-01) to Q1: which settings or input produced the 25-fly per0 DD reference set?
Answered during the run, relayed through Sheyum; recorded here for the file.
- Input: not the raw *_from1stLD.txt files but the pre-cut copies *_from1stLD.txtcutdata_4.8.2025.txt
  (per0DD_activity_files_1sec_res.zip, SHA-256 prefix 778599db; 16 x 4.000 d + 9 x 3.000 d = 91.0 fly-days;
  the 9 short ones are exactly the runs 2_18_2018 and 6_25_2018_2boardmb). Sheyum: the cut removed the first
  24 h of DD. The CS DD files were cut the same way before export.
- Settings: never recorded. With the cut files, Teliminate 0 gives the right 25 flies; neither the current
  nor the pre-09-22 pipeline reproduces 3989 / 3987 (both give 4017 / 4022, bit-identical; edge bouts dropped).
  The stored file is lost, so its exact provenance is unresolved.
- Resolution: the replication route. The re-extracted set per0DD_bouts_20r_replication.mat (SHA-256 b14ba340...)
  is the per0 DD reference from now on; the 2026-09-15 rise test replicated on it (root ledger, 20r entries).
