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
