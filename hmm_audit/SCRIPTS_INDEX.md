# Script index — for converting the analysis to MATLAB

**For the S/Q/W analysis (PRESPEC 19-20e), port `sqw_ref.py` following `PORTING_GUIDE.md`, tested against `golden/`.** The table below is provenance.

Python 3.13; packages: numpy, scipy, matplotlib, hmmlearn 0.3.3, scikit-learn.
Every script hard-codes the session scratchpad path
`/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/...` for its data folder
(`CD`, `D` or `sys.path` lines near the top). Change those to local paths before running.

Data layout expected by the count-data scripts: `counts/<group>/**/*.txt`, one integer per line =
movements per second, starting at ZT0, groups `CS`, `fumin`, `per0LD`, `per0DD`, `sss`.
The scored-data scripts read the `.mat` cell arrays (`C`: model index, P 3x3, B 3x2 per fly-window; `R`:
observation sequences, 1 = inactive, 2 = active) through `common.py`.

Plan files (`PRESPEC_*.md`) give each analysis's definitions. `RESULTS_LEDGER.md` gives every result,
including failures. Use them as the specification when porting.

## Current analyses behind the figures (priority for conversion)

| Script | Plan | What it computes | Outputs |
|---|---|---|---|
| `e23.py` | 15, 16 | E2: activity after waking vs before sleep (min 2-10), bout-swap permutation test. E3: run-triple time-irreversibility (KL vs time reversal), segment-reversal null. Args: group names. | `e23.log`, `e23_new.log` |
| `taper.py` | 17 | Per-fly asymmetry D, sign tests, taper depth; part C (relaxation fits) is INVALID, kept for the record. Defines `fly_bouts`. | `taper.json`, `taper.log` |
| `figs.py` | 17 | Figures 1-4 (2-10 min window, pre-registered); Fig. 3 exit depth and Fig. 4 are descriptive. Imports `taper.py`. | `figures/fig1-4_*.pdf/png` |
| `v18.py` | 18 | Sleep scoring with 1 s bridging (lab rule: only between two >= 300 s inactive runs); versions V0, V1b, V2b (+-5 min). Defines `sleep_runs(x, bridge)`. | `v18b.json`, `v18b.log` |
| `figs_pm5.py` | 18 | Figures 1-3 for the +-5 min window (descriptive), bridged scoring. | `figures/fig*_pm5_bridged.*` |
| `cycle.py` | 19 | Signed cycle current: Q/A classification at each transition (theta = 0.33/0.5/0.75 x fly mid-wake activity), c = P(Q|onset) - P(Q|waking), J per hour, bootstrap CIs, sign/Wilcoxon/Mann-Whitney. Imports `v18.py`. | `cycle.json`, `cycle.log` |
| `fig_cycle.py` | 19 | Figure 5 (cycle diagram, c, J). | `figures/fig5_cycle_current.*` |
| `epr_min.py` | 15 | Minimum and maximum EPR and cycle flux over the equal-likelihood family of the scored 3-state model. | `epr_min.json`, `epr_min.log` |

## S/Q/W HMM (PRESPEC 20-20e)

| Script | Plan | What it computes |
|---|---|---|
| `sqw_ref.py` | 19-20e | REFERENCE implementation, one function per step (port this) |
| `verify_ref.py`, `make_golden.py` | -- | checks sqw_ref against the pipeline; writes `golden/` test vectors |
| `sqw.py`, `sqw_common.py`, `sqw_synth.py`, `sqw_power.py` | 20 | per-window design (failed V1), steps-per-fit check |
| `sqw2_common.py`, `sqw2.py`, `sqw2_synth.py`, `sqw2_report.py`, `sqw2_countJ.py` | 20b | per fly x light/dark fits, routing, count J, Jensen, rule-of-3 bound |
| `sqw3.py` | 20c | 6 h blocks |
| `sqw4.py` | 20d | reversible null, attenuation |
| `sqw5.py` | 20e | across-day and within-day direction flips |

## Scored-HMM audit (identifiability, decoders)

| Script | Plan | What it computes |
|---|---|---|
| `common.py` | — | Loads `.mat` files; accepted windows; Xn_pred / Yn_pred formulas |
| `analysis_AB.py` | 1 | Emission degeneracy; identified set of Xn_pred, Yn_pred over S P S^-1 |
| `analysis_C.py` | 1 | Viterbi vs FFBS route counts per window (FFBS in Python) |
| `analysis_D.py` | 2 | Joint identified set of (Xn, Yn), Delta, Sigma |
| `plot_scatter.py` | 4 | Per-fly Xn vs Yn scatter with identified ranges |
| `viterbi_vs_ffbs.py`, `equivalent_models.py` | — | Synthetic demonstrations (no fly data) |

## Bout-level and count-HMM analyses (supporting, mostly superseded for the paper)

`pooled.py`, `pooled_LD.py` (PRESPEC 5/5b), `hetero.py` (6), `within_fly.py`, `within_fly9.py`, `em_fit.py`
(7, 9), `irrev.py` (8, 8b), `counts_analysis.py`, `counts_C3_resume.py` (10), `flatK.py` (11),
`hhmm.py`, `hhmm_day.py`, `hhmm_fly.py`, `synth_test.py`, `validate13.py`, `validate13_day.py` (12-14),
`aging_check.py` (descriptive).

## Known pitfalls to carry into the MATLAB versions
- Exclude minute 1 next to every transition (a movement is guaranteed there by definition).
- Never normalise a bout by its own mid-wake activity (near-zero divisors; PRESPEC 17 C).
- Cut time windows at fixed clock times, never extend them to finish a bout (length bias; PRESPEC 8b).
- The +-5 min window is post hoc; the 2-10 min window is the pre-registered one.
