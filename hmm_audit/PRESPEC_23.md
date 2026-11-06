# PRESPEC 23 — Shaker minisleep (Sh^mns): which way does a third short sleeper run the cycle?

Date 2026-11-06. Written and committed BEFORE any behavioural quantity of these flies has been computed. Only file
metadata (count, names, lengths) has been inspected. This is an exception to the freeze (new data), approved by Sheyum.

## Data
shk-mns_LD.zip (SHA-256 9697ee5a...): 64 files, male Sh^mns, backcrossed to CS (as fumin and sss were), LD 12:12, files
start at ZT0. They come from several recording runs (2018 dates in the names); lengths are 2.5-9.3 days. Raw data are
not uploaded to the repository.

## Records and readouts (identical to the other genotypes)
- Death rule (Sheyum's, as PRESPEC 24): the first run of >= 12 h with no activity marks death; the record is cut 24 h
  before its start.
- Days: up to the first 4 whole days (the CS / per0 record length); at least 2 whole days, else the fly is excluded.
- Symbols sqw_symbols (theta 0.5, light/dark references), 12-h phases, mc_fit chain, observable_quantities (c_obs,
  count J, count EPR, chain J / A / EPR with the rule of three, landscape), sleep min/h.
- No-cycle fly-phases (|J_chain| < 1e-12) are counted and excluded from the direction tests. Values are rounded to
  10 decimals before ranking.

## Why: the act 3 question
sss reverses the cycle by day; fumin runs it forward, hardest of all; sleep amount does not set the direction (C14).
SLEEPLESS regulates Shaker. Sh^mns asks whether the sss reversal follows the Shaker pathway.

## PRIMARY (Holm over 2)
P1. Direction by day: Wilcoxon signed-rank of J_chain against 0 over light-phase fly-phases with a cycle, two-sided.
    Reported with the fraction of flies with J < 0.
P2. Direction at night: the same test in the dark phase.

## SECONDARY
S1. C14 within Sh^mns (Holm over 4): Spearman of sleep with c_obs and with J (cycle flies), light and dark, 10,000
    permutations. The sign of the J-sleep slope by day says which way extra sleep turns the cycle.
S2. Onsager conductance L (EPR = J^2/L fit, as C5) against CS: fly-bootstrap difference, two-sided.
S3. Landscape: the sign of E_S - E_Q, by day and night (median with fly-bootstrap CI).
DESCRIPTIVE: matched-sleep comparison with fumin and sss (the C14 overlap method); the cumulative walk; HMM level only
if Sheyum asks (post hoc).

## Reading (fixed now)
- "Sh^mns reverses the cycle by day, like sss": P1 Holm < 0.05 with median J < 0.
- "Sh^mns runs forward by day, like fumin": P1 Holm < 0.05 with median J > 0.
- Otherwise: "no net direction by day detected". Then the fraction negative and the CI are reported, never "intermediate".
- Comparisons with CS, fumin and sss cross recording runs and years (as for fumin and sss): stated with each.
- No other day count, phase split or readout is added after looking.

## Output
`sqw23.py` -> `sqw_out/s23_shmns.json`, `figures/fig23_shmns.png`. Ledger entry with commit hashes.
