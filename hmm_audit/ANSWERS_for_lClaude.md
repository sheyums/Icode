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

## A2 (2026-11-01) to Q2: stored 20b fits
Committed under hmm_audit/: sqw20b_fit2_theta0.5.jsonl (= sqw_out/fit2_0.5.jsonl: per fly g, fly index in sorted-glob file
order, T, ref, frac, F{K: ll, pi, P (nP x K x K), B, it} for K = 2, 3, 4, and the labelled K = 3 fit P3, B3, pi3, deg, xn, yn)
and sqw20b_route2_theta0.5.json (= sqw_out/route2_0.5.json: key "g|fly"; ffbs/vit/pairs/hours, per-path pair counts
'pairs_path' (200 paths x phase x 3 x 3) and route counts 'routes_path'). All indices 0-based; Q/W/S = 0/1/2.
Note for the 1e-6 check: the stored ll is of the BEST of 20 numpy-seeded starts; reproduce it by evaluating hmmForward at the
stored (pi, P, B) (exact), and by your own 20-start fit reaching it within tolerance (statistical), not by matching starts.

## Answers to the five cases in PORT_VALIDATION.md stage 1
1. sss: the paper's set is all 57 files. The three you found unused are sss flies 34, 35, 50 (4_20_2018 #5, 5_9_2018 #3,
   9_16_2019 rc2 #3): the three flies with NO Q/W reference in either phase (all wake minutes missing). They are excluded by
   rule from everything needing Q/W (20i; 20p amendment B), and are also the three HMM-degenerate sss fits. Nothing to change.
2. Agreed: "no cycle" lumped three kinds. Correction appended to RESULTS_LEDGER.md (20s-A correction); the Methods wording
   will name the kinds, as 20q amendment E does.
3. Agreed: dropping Inf biases the chain-EPR median low exactly where irreversibility is largest. Correction appended:
   chain-EPR medians are reported with Inf ranked highest (your EPR_chain_med_withInf), with the Inf share beside them.
4. Fine as implemented.
5. Thank you: the calibration (Wald 22% false positives at nominal 5%; bootstrap difference 5.2%) is now cited in the
   ledger as support for using the bootstrap difference.
