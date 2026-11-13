# beam_hmm — sleep-wake dynamics from 0.2-s multibeam position

New project, separate from hmm_audit (whose S/Q/W paper is under an analysis freeze since 2026-11-01).
Root CLAUDE.md rules 1-5 apply here: PRESPEC committed before computing, append-only ledger with commit
hashes, no re-binning or widening after looking, narrowest block plus a sample-size-matched null for
every reported number.

- `PRESPEC.md` — one entry per analysis, committed before it runs.
- `RESULTS_LEDGER.md` — every result, including nulls; append-only.
- `data/` — raw files, local only (gitignored). Never commit raw data.

Prior history that shapes the design: hmm_audit PRESPEC 12-14 and DESIGN_REVIEW.md (hierarchical count
HMMs failed synthetic validation: per-fly routing not estimable; pooling made a false state from fly
heterogeneity; free per-fly emissions made states degenerate).

## Tools available elsewhere (not in this repo)
- `EstimatePhaseStamp.m` (user's MATLAB analysis folder; tested there with other scripts): recovers startZT per
  experiment from pooled activity. One anchor peak, partner derived from the photoperiod, lights-off labelled
  by a 4-measure vote (ClockState "intact") or by the larger rise ("arrhythmic"). 10-min bins for intact flies,
  first cycle skipped. Candidate uses here: a per-file check that row 1 really is ZT0 before phase-dependent
  analyses (a wrong phase is a silent 12-h inversion); for DD data, the anchor idea combined with a period
  estimate (e.g. chi2p), since folding at 24 h fails when the clock free-runs. Not needed while the files'
  stated start (lights-on) is trusted. Note: its arguments-block comment for BinMinutes says 15/60, the code
  and help say 10/30.
