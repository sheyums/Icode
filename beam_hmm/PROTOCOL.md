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
