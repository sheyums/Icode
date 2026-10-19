# HMM routing audit (separate from the Icode bout-fitting suite)

Audit of `compute_HMM_energetics_test.m` / FFBS route counts on the CS, per0, fumin and sss
3-state HMM fits. Not part of the closed bout-duration library; nothing here is called by it.
The raw `.mat` data are NOT included.

- `PRESPEC*.md`: analysis plans, each committed before its analysis ran.
  `AUDIT_GIT_LOG.txt` is the commit log of the local audit repository that carried those
  timestamps (it lived in an ephemeral scratchpad and could not be pushed directly).
- `RESULTS_LEDGER.md`: append-only results, including crashes, nulls and abandoned analyses.
- `common.py`, `analysis_*.py`, `pooled*.py`, `hetero.py`, `within_fly.py`: Python
  reimplementations (no MATLAB was available); `*.log`: their outputs.
- `viterbi_vs_ffbs.py`, `equivalent_models.py`: synthetic demonstrations (no fly data).
- Paths in the scripts point at the original session scratchpad and must be edited to rerun.
