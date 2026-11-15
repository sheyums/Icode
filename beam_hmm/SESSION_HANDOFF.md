# beam_hmm — session handoff (written 2026-11-15, before context compaction)

Read this first after a compaction, then PLAN.md, PRESPEC.md (entries), RESULTS_LEDGER.md (latest entries at the end).
Branch: `claude/funny-johnson-5t075n` (repo sheyums/Icode). Everything below is committed and pushed.

## Who / what
- User: academic physicist/biologist; wants answers grounded in published work. Rules in repo-root CLAUDE.md apply
  (PRESPEC committed BEFORE computing; append-only ledger with commit hashes; no re-binning after looking; narrowest
  block + matched null; flies are the unit). Disclose anything decided after seeing data; label post hoc.
- Project: richer model of Drosophila sleep-wake from multibeam position data (brief in the first user message:
  hierarchical HMM vs time-dependent HMM; irreversibility, homeostasis, circadian decline, approach to death).
- Data: ONE pilot run, CS males, LD 12:12, 9 weekly raw files (beam_hmm/data/week1..9, gitignored; .npz made by
  explore/load_week.py). Only channels 17-32 readable. A larger repeat experiment is planned (months away).
  Format: 3 header lines; 64 positions + time; every other line is 64 x -1 (dropped); true step 0.4 s.
  Position 0 = "no beam broken" (beyond beam 1 in live flies); beam 15 dead (14<->16 reads as a 2-beam step).

## Pipeline facts (fixed by PRESPECs P1, P1-B/C/D/E/F)
- Spike filter (A-B-A, |B-A| >= 4 -> A). Dead-hub repair (hub = top target of >= 4-beam jumps with share >= 30% and
  returns >= 20%; every hub reading replaced by the last other reading). Narrow-range rule (< 8 beams spanned) -> exclude.
- ZT from the flies' own lights-off rise per file (P1-E): files 1, 2, 5, 7 start 30-40 min from lights-on (clock drift).
  Week 6 by hand check: lights-on 23.83 h into the file. First 24 h of each file excluded; transfers = gaps between files;
  week 1 from hour 150 excluded (tubes out ~153 h).
- Deaths: P1-F (DEAD-LIKE < 100 changes/week; UNREADABLE = hub > 50% or narrow; outages are not deaths). Matches the
  survival sheet (fly No. = channel number) except ch25 (undetermined). ch27, ch30 unreadable (alive per sheet).
- 14 fly-weeks: live flies HIDDEN behind a blocked beam 16 (constant "16" all week) — plug/food on the end beam.
- Lab sleep scoring (sleeptrace + ironout) == beam_hmm scoring second for second (0 of 28M s differ).

## Results so far (all in RESULTS_LEDGER.md with hashes)
Real data (pilot; weeks 2-3 = unseen test set; week 1 descriptive):
- A1 SUPPORTED: by day, sleep at tube ends more than short pauses (11/12, Holm 0.0067); not at night.
- A4 SUPPORTED: slower walking at night (0/10 positive, Holm 0.015).
- A2 not supported (pauses before sleep); A3 not supported (flicker day/night) -> flicker likely detector effect;
  sleep measure robust to flicker.
- A2b (port of hmm_audit PRESPEC 21): primary verdict "no shape asymmetry"; secondary 5-s setting strongly asymmetric,
  OPPOSITE direction to the earlier project (less active first ~25 s after waking). Ruled out: port error (exact
  reproduction of per0 DD, repro_a2b_per0dd.py), detector, lab conversion, scoring, 1-h switch exclusion.
  Open: cohort/age/genotype/LD-vs-DD (PLAN.md "Open question").
- A5 ageing: no trend survives Holm (sleep, end-resting, walking speed nominal, same direction).
- T1-C1 SUPPORTED (2026-11-14): shared time-of-day modulation of observed sleep onsets and wake-ups beyond the L/D step
  (11 flies, LR 1012, rotation p 0.001; robust to exclusion 1/2/3 h). Onsets peak ZT 5.1 (siesta) and 16.9; wake-ups
  peak ZT 10.8 (evening anticipation) and near lights-off. LD = light + clock.
Synthetic gates:
- D1-B FAIL (windows -> transitional states); D1-D FAIL (degeneracy rule forced K = 2; serial correlation makes spurious
  states as large as real ones); T1 FAIL (per-fly clock power); T1-B PASS (cohort clock test, hidden states, K fixed);
  T1-C FAIL (3/20 false); T1-C1 PASS (2/60 false, 20/20 power).
- B1 (direction-aware AR-HMM at 0.4 s, = D1-E) gate RUNNING on the user's machine (see below); breadth run (7 days,
  10 flies) queued after it.

## Running elsewhere: the user's local Claude session
- Session "Set up Python and run beam_hmm synthetic gates", id session_01LXkmseoPDPi3ntdCxr7Ehn (Windows, Remote
  Control). I can read its transcript (mcp list_events) and send_message to it; it CANNOT message back (it posts in its
  own transcript). Scheduled check-ins (send_later) from this session fail with access denied -- don't rely on them.
- Clone: C:/Users/sheyu/beam_hmm_work/Icode; venv python at ../../venv/Scripts/python.exe; raw files at
  C:\Users\sheyu\Box\Data from Zichao\Super long multibeam exp\raw_positions_all64_weekN.txt; MATLAB R2026a available
  (headless -batch runs hung on exit once; LoadingCheck now closes saved invisible figures).
- LAST REQUEST SENT (2026-11-15): restart the B1 gate with OMP/OPENBLAS/MKL/NUMEXPR threads = 1 (thread over-
  subscription slowed gates ~50x here); report the speed after ~15 min; push results/b1_gate + REPORT.txt when done.
  When B1 reports: log it in the ledger (verdict from the files), then follow the plan below.

## MATLAB tools for the lab (beam_hmm/matlab; all with TYPICAL SEQUENCE in the header)
LoadingCheck (2-h re-seat screen: BLOCKED / FROZEN / NARROW / LOW ACTIVITY, per-channel figure, WHAT TO DO printout;
verified in MATLAB, figure in results/loadingcheck), BeamFingerprint (QC + 4 figures), RawToActivity1s (raw -> beams per
second, 'beam_hmm' or 'legacy' mode, QC, refuses to append files), golden tests test_BeamFingerprint /
test_RawToActivity1s (full-week MATLAB runs still pending). Python twins in labconv.py. Octave 8.4 is installed here for
logic checks (no graphics, no `contains`/`gobjects`/`table`).

## Next steps (decided / proposed)
1. When B1's gate reports: PASS -> T1-B on real flies with B1 states; then B2 (fast hidden states nested inside the
   OBSERVED S/Q/W classes). FAIL -> no likelihood-based state discovery at 0.4 s; options: fixed K + stability, semi-
   Markov models, bout-duration suite; T1-C1 stays the clock result.
2. Offered, not yet chosen: T1-C1 descriptively across weeks 4-9 (circadian control vs age, for A5); design note for
   the repeat experiment (flies per question from A-measure spreads; LoadingCheck; record food end per tube; DD segment).
3. Open questions to the user: monitor manual (what is reported when two beams are broken; meaning of 0); food end per
   tube for this run (to test food-end vs plug-end blockages).

## Working habits that mattered
- Run multiprocessing gates with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1.
- Keep this session active while computing here (idle containers are reclaimed; background jobs die).
- Commit code before running; commit results after; append-only ledger, corrections appended.
- Use `git pull --rebase` before pushing (the local session pushes results to the same branch).

## Update 2026-11-15 (after compaction)
- B1 main gate on the user's PC: ~15-18 core-h per fly (not 3-4 h as I estimated here); 3/40 done at 12:56 local;
  finish expected ~midday 11-17 or later. Breadth run DROPPED by the user (no b1_gate7). Local session pushes
  results/b1_gate + REPORT.txt when the launcher exits.
- DESIGN_NOTE.md written (repeat experiment: ~70 loaded flies recommended; QC protocol; DD segment; compute budget).
