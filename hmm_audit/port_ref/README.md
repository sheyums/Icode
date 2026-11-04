# Reference outputs for the MATLAB port (derived results only; no raw data)
Copied from hmm_audit/sqw_out/ (which is not synced). Python versions as of 2026-11-03.
| file | produced by | what to compare |
|---|---|---|
| s22_sleep_vs_J.json | sqw22.py (PRESPEC 22) | per fly-phase sleep, c_obs, J; per-block rho (exact), p (MC, approximate), Holm |
| s22h_hidden.json | sqw22h.py (post hoc) | hidden sleep (60*pi_S of P3), J_hmm, d; per-block rho, p, Holm |
| s25_windows.json | sqw25.py (PRESPEC 25) | per-window N_SW, N_WS, class (exact); null means (MC, approximate) |
| s26_loopseq.json | sqw26.py (PRESPEC 26) | per-fly n_loops, S_obs, run lengths (exact); S_null, E, p (MC, approximate) |
| cycle_times.json | cycle_times.py (descriptive) | loops/h fwd/back, h per net cycle, loop durations (exact) |
| r3_model.json | sqw_r3.py (rule-of-three amendment, revised 663a0c2) | A_r3, EPR_r3, landscape, HMM model EPR per fly-phase |
2026-11-04: s22/s22h re-generated with ranks of round(x, 10) (ties that differ only by round-off now tie).
