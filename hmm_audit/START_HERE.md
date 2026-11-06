# START HERE

**Reading to understand the science (added 2026-11-06), not to implement it:**
1. The report, "The Measured Sleep-Wake Cycle" (claude.ai artifact 2JFn62CYxDRbgfn8pTzpcZ): every finding in one
   sentence and one figure, with a physics and a neurobiology reading. Its first section, "The story", is the paper's
   three-act structure.
2. `PLANS.md`: "The story" (the acts, with the HMM's role in each) and the claim list C1-C14.
3. `DRAFT_HMM_motivation.md`, `DRAFT_act2_loops.md` and `DRAFT_act3_C14.md`: draft paper text.
4. `RESULTS_LEDGER.md`: the record behind every number. Look up entries by PRESPEC number; it is long and
   append-only, so later corrections follow earlier entries.
Time-reversal asymmetry in brief: C2-C5 (the directed W->Q->S->W cycle, flux, affinity, entropy production), C9/C9*
(irreversibility beneath the observed chain, dwell times that depend on the direction of travel) and C13 (the
asymmetry in raw 1-s activity: flies wind down gradually and wake abruptly).

---

# Implementation guide (for lClaude or anyone writing the MATLAB functions) — state as of 2026-10-31

This folder holds several weeks of analysis history. Most of it is RECORD, not instruction. For the paper (Plan A)
read only the files marked CURRENT, in this order. Where an older file disagrees with a CURRENT one, the CURRENT one
wins.

## 1. What the paper claims and why
| file | status | use |
|---|---|---|
| `PLANS.md` | CURRENT | ONE merged paper (2026-11-01): claim list C1-C11, each with the analysis that supports it, and the paper outline. Plan B group-level tests 1-4 are on branch planB-feasibility (RESULTS_LEDGER_planB.md). |
| `RESULTS_LEDGER.md` | CURRENT (append-only record) | every number, with the commit that produced it. Read the entries named in PLANS.md, not the whole file. Earlier entries are history; corrections are appended, never edited in place. |

## 2. What to implement in MATLAB
| file | status | use |
|---|---|---|
| `sqw_ref.py` | CURRENT, the ONLY code to port | one function per step. Sections 1-14: symbols, HMM, counting, round trips, dwell, Onsager, route-vs-duration, landscape, 30-s grid. **Section 15: Plan A** (observed chain, block irreversibility, TUR, observable quantities). |
| `golden/G1..G12.json` | CURRENT | fixed inputs and expected outputs. Port each function until its golden test passes. G12 covers section 15. |
| `PORTING_GUIDE.md` | CURRENT | function map, conventions, algorithm details; **section 11 = Plan A** (read it before section 15 of sqw_ref). |
| `FIGURE_RECIPES.md` | CURRENT | data extraction and plotting per figure; **section K = the Plan A figures A1-A4**. Sections A-J remain valid for the supplementary figures. |
| `verify_ref.py`, `verify_planA.py` | CURRENT (Python only) | show that sqw_ref reproduces the stored results (58/58 for Plan A). They need the raw data; not to be ported. |

The MATLAB side already has `epr_count_based.m` (identical to `epr_count_bound`; reuse it) and the bout suite in the
Icode root (`BoutHazard`, `CompareBoutModels`, ...), governed by the root `CLAUDE.md`.

## 3. Definitions of each analysis (the plans)
CURRENT for Plan A: `PRESPEC_20.md` (+ amendment A), `20b`, `20d`, `20g`, `20h`, `20k`, `20m`, `20n`, `20o`, `20p`
(+ amendments A, B), `20q` (+ amendments C, D, E), `20s.md` (Plan A observable recomputation). Supporting, still valid:
`20e`, `20f`, `20i`. Superseded for the paper's headline numbers but kept as the HMM bracket: `20c`, `20j`, `20l`
(their observable versions are in `20s`).
The repo-root `PRESPEC.md` holds the MATLAB bout-suite plans, including **20r** (hidden steps, CS DD; instructions in
`INSTRUCTIONS_20r_lClaude.md`).

## 4. HISTORY — do not port, do not follow
- `PRESPEC_2` ... `PRESPEC_19` and their scripts: the old scored S/W per-window HMM audit, count HMMs, hierarchical
  HMMs (PRESPEC 12-14 failed validation), bout-level analyses. They explain how we got here.
- `PRESPEC_20.md`'s per-window design (V1) failed validation; the per-fly design of `20b` replaced it.
- Analysis scripts `sqw*.py`, `fig_*.py`, `*_report.py`: they contain paths, multiprocessing and the history of each
  run. Their logic for the paper is in `sqw_ref.py`; use the scripts only to see how a figure was assembled.
- `SCRIPTS_INDEX.md` lists every script with its status.

## 5. Two things that are easy to get wrong (details in PORTING_GUIDE section 11)
- The TUR bound is **B = J^2 / D** with D = Var(N_t)/(2t). (A plan once had 2 J^2/D; corrected before any data run.)
- Every comparison of observable vs hidden quantities is made **per fly x phase**, then summarised over flies.
