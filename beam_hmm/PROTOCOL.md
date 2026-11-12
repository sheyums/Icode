# Protocol facts (as given by the user; what is still unknown is listed as unknown)

Source: user's info sheet "Super long CS multibeam exp" (PDF, uploaded 2026-11-12), plus chat.

## Stated
- Genotype: CS (Canton-S). Multibeam tubes, 16 beams per fly; one column per fly, one row per ~0.2 s sample.
- 9 raw files, each roughly one week. Every file begins at lights-on (row 1 = ZT0 of that file).
- Each fly keeps the same column position in every file.
- Food end differs between flies: near beam 1 for some, near beam 16 for others. Which end, per fly, is
  NOT recorded. The user's prior: flies typically spend the midday nap near the food end.
- Survival as recorded (coarse, by week):
  - No. 1 died in/after week 4.
  - No. 10, 14, 25, 41, 47, 50, 58 died in/after week 5.
  - No. 2, 4, 8, 12, 13, 15, 16, 19, 20, 23, 24, 27, 30, 31, 32, 33, 34, 36, 38, 43, 44, 48, 53, 57, 61
    alive in/after week 6 (25 flies).
  - No. 4, 12, 19 and one more fly (position unknown) alive in/after week 8.
- Flies were 3-4 days old at the start (brief).

## Unknown -- to be confirmed or measured from the data
- Total number of flies loaded and the position numbering (the highest listed is 61).
- Light schedule (12:12 LD?) and whether it was constant over 9 weeks; temperature; sex.
- Whether the 9 files are contiguous in time or separated by gaps (food changes, tube transfers).
  A transfer or restart is a break in the record.
- Exact death times: to be inferred from the data (silence >= 12 h; cut 24 h before), with the
  sheet's weekly list used only as a consistency check.

## Consequence fixed now (before seeing any data)
Food end must be inferred per fly. It must NOT be inferred from midday-sleep position and then used to
report that flies nap near food: that is circular. Any inference rule is pre-specified in PRESPEC.md,
uses a block that is then excluded from analyses involving place, or uses a criterion independent of
the quantity being reported.

## Added 2026-11-12 (user, chat)
- All flies male. LD 12:12 throughout.
- Roughly every 7 days flies were moved to a tube with fresh food and put back in the monitor: every file
  boundary is a transfer (a disturbance and a break in the record), not a continuation.
- Death rule (lab convention): no activity for >= 12 h = dead; exclude the preceding 24 h as well.

## File format (week 1, measured on the file itself)
- 3 header lines: "Recording started at: Wed Sep 09 04:26:13 2026", blank, "Recording frequency: 5 Hz".
  CRLF line endings, tab-separated.
- 3,024,000 data lines = 604,800 s x 5 Hz, 65 fields each: 64 channels then a time column.
- Lines alternate: a line of beam positions with a time stamp, then a line of 64 x "-1" with an empty time
  field. Meaning of the -1 lines is to be asked, not assumed.

## Lab preprocessing (PreProcess.m, supplied by the user 2026-11-12) and what week 1 says about it
PreProcess reads every line after the header (including the all -1 lines) and, per column, replaces
-1 and 0 by the previous value, and 1 by the previous value unless the previous value is 2. Leading 0/1 are
set to the first other value.
Measured on week 1 (channels 17-32, qc5.py):
- Every run of 0 is entered from 1 (100% in every channel) and mostly left to 1 or 2. So 0 is a POSITION,
  "past beam 1" (the tube end beyond the first beam), not a dropout. PreProcess folds it into beam 1.
  Share of time at 0: up to 9.4% (ch17); 6-7% in ch22 and ch24, which sit at low beams by day.
- No equivalent code exists at the beam-16 end.
- Value 1 is mostly entered from 0 or 2 (legitimate); entries from 3+ are rejected by PreProcess, which acts
  as a jump filter at that end only. Hub artifacts at other values (3, 5, 11, 16, ...) are not touched.
- Because the -1 lines are kept and filled, PreProcess output has every 0.4-s reading duplicated at 0.2 s.
- PreProcess does not explain why channels 1-16 read only 13-16.
Consequences for this project (fixed before any analysis): keep 0 as its own position state ("beyond
beam 1") rather than merging it; drop the -1 lines rather than fill them (true step 0.4 s); any speed is
computed per 0.4 s. Detection is place-dependent: movement beyond beam 1 is unseen, so stillness there is an
upper bound.
