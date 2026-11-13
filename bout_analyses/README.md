# bout_analyses: MATLAB scripts behind the 2026-11 ledger entries

These scripts were written and run in a Claude Code session's scratch folder
between 2026-11-04 and 2026-11-12; the ledger entries call them "session
scratch, not in the repo". They are committed here unchanged so every entry
points at code that exists. They are records, not a maintained library:

* Paths inside them are absolute paths to that scratch folder
  (`C:\Users\sheyu\AppData\Local\Temp\claude\...\scratchpad\...`), including the
  staged data folders (`per0dd_stage`, `csld_stage`, `csdd_stage`,
  `tim0dd_stage`) and the intermediate `.mat` files, which are NOT in the repo.
  To rerun one, point those paths at your own copies.
* Data inputs: the bout tables and raw activity files are in the Python
  reproduction kit (`bout_reproduction_kit_2026-11-12.zip`, analysis folder),
  which also holds the MATLAB reference results and a Python port.
* Scoring used the analysis-folder copies of `process_LocomToSleep` and
  `extractBoutsForModelCompetition` (analysis repo commit 04e0508).

| folder / file | ledger entry | what |
|---|---|---|
| `homeo_run/` (`run_homeo`, `run_B`, `run_E`) | homeostat entries A-E (2026-11-06) | builds `homeo_results.mat` (all bout sets); cross-type coupling, cumulative-sleep tests |
| `per0dd_run/`, `csld_run/` | time-arrow entries (2026-11-05) | BoutIrreversibility on per0 DD and CS LD |
| `ar1_run/` | AR(1) noise entry (2026-11-06) | BoutAR1Noise / residual families |
| `skew/` | GH noise-skew checks (2026-11-04/05) | |
| `F_run/` (`run_F`, `run_G`, `run_H`, `fig_H`) | F, G, H (2026-11-08/09) | Shreenath's exports; extreme-bout clustering vs his fitted replicates |
| `I_run/`, `J_run/`, `K_run/`, `L_run/` | I, J, K, L | clustering source ladder; day x ZT null; trend vs fluctuation; per0 after day 1 |
| `M_run/` (`run_M`, `diag_M`, `fig_M`, `fig_M2`) | M (+ correction 2026-11-10) | sleep/wake neighbour coupling, long-wake tail |
| `recal/` (`calK`, `calL`, `calM1`, `calM2`) | correction 2026-11-10 | calibrations rerun with independent seeds; M counter fixed (`calM` is the superseded first attempt) |
| `N_run/` | N | tim0 DD first-day check (also builds `tim0_flies.mat`) |
| `O_run/` | O | tim0 DD repeats of L and M |
| `P_run/` | P | CS DD coupling |
| `Q_run/` (`q*.m`, `run_Qcal`, `run_Qdata`) | Q | one-step link, kappa |
| `bout_kit_build/export_kit.m` | - | exports the Python kit's data and reference values |
| `stage.m`, `stage_cs.m`, `dupcheck.m`, `inv2.m`, `homeo_run/stage_check.m` | - | data staging and duplicate checks |
| `binomcheck.m`, `bc/`, `diag_strata.m`, `eps_per0ld.m`, `frf_try.m`, `ghtime.m`, `per0ld_phase.m`, `phase_cs.m`, `I_run/occ.m` | - | one-off diagnostics referred to in the ledger notes |
