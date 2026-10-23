# Reminders (open items)

- [ ] **Threshold-free S/Q/W HMM** (user request, 2026-10-23): after the thresholded S/Q/W thread is complete,
      fit a version without the Q/W threshold: each wake minute carries its number of active seconds (0-60) and
      the HMM learns the wake states from binomial-type emissions, so H1 emerges as "low-activity wake" with no
      theta. Compare its routing with the thresholded version.
- [ ] Convert the Python scripts behind the figures to MATLAB functions (see SCRIPTS_INDEX.md).
- [ ] Supplementary versions of Figs 2 and 4 with all genotype x light/dark blocks.
- [ ] E3 irreversibility as a per-hour lower bound (careful unit conversion).
- [ ] Independent sleepless flies to test the reversed cycle current (PRESPEC 19 is descriptive).
