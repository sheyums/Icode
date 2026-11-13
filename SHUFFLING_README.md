# Shuffling detection

Detects sustained alternation between two adjacent beams (`15,14,15,14,...`)
in the cleaned position trace, at native 5 Hz resolution.

The point of working here rather than downstream is that
`multibeam_locomotion` computes `sum(abs(diff(position)))`, which conflates
path length with net displacement by construction — oscillating in place and
walking steadily give the same number. That distinction only exists before
the collapse to a locomotor scalar.

The same detector does double duty as QC. Near-continuous alternation across
a whole recording is a hardware fault, and the existing death detection
structurally cannot catch it: that looks for *inactivity*, and this is
continuous *activity*.

## Files

| File | Role |
|---|---|
| `find_shuffles.m` | Core detector. One position vector in, episodes out. |
| `PreProcessV2.m` | Non-interactive, crash-guarded copy of `PreProcess.m`. |
| `shuffling_qc.m` | Per-fly status table (stage 1). |
| `detect_shuffling.m` | Driver: clean → QC → detect → chattering → sleep overlap. |
| `shuffling_stats.m` | Population views, per alternation-length threshold. |
| `validate_shuffling.m` | The two post-build checks. |
| `test_shuffling_detection.m` | 22 tests. Run from this directory. |

Nothing here modifies `PreProcess.m`, `multibeam_locomotion.m`,
`process_LocomToSleep.m`, `sleeptrace.m` or `segment_sleep_phases.m`.

## Quick start

```matlab
r  = detect_shuffling('rawfile.txt', struct('N', 64, 'Verbose', true));
st = shuffling_stats(r, struct('Thresholds', [3 5 10]));
```

Add `'OutputDir', 'out'` to write the QC and episode tables as
tab-delimited `.txt`.

## What still needs wiring up

Three things could not be built here because the functions they depend on
were not available. Each is a **declared input**, not a stub — the code
paths that consume them are written and tested against synthetic data, and
they will start producing real numbers as soon as they are fed.

### 1. The analysis window — the highest-risk input

`process_LocomToSleep` trims (`Teliminate`) and truncates (death detection)
before computing anything, so `sleep_out` indices are **not** raw file rows.
Pass the same window it used:

```matlab
opts.Window = [firstSample lastSample];   % 64 x 2, raw file row indices
```

Omit it and whole files are analysed, every fly is flagged
`WindowNotSupplied`, and a warning fires. Nothing errors — that is exactly
the silent misalignment to avoid. Prefer having `process_LocomToSleep` hand
the window over directly to re-deriving it.

Every episode carries **both** `StartSample` (raw file row) and
`StartSampleInWindow`, so the two frames can never be confused by accident.

### 2. The sleep/wake trace

```matlab
opts.Sleep = struct( ...
    'Trace',       sleepOut, ...      % bins x 64, or 1x64 cell of vectors
    'BinSeconds',  60, ...
    'StartSample', firstSample, ...   % raw row that bin 1 corresponds to
    'Polarity',    'ZeroIsSleep');    % REQUIRED
```

`Polarity` deliberately has **no default**. `sleeptrace` returns an activity
wave in which `0` means sleep — the opposite of what the name suggests, and
the gotcha its own docstring and `process_LocomToSleep`'s both warn about. A
wrong guess would invert every sleep statistic silently, so it has to be
stated. Use the cell form for `Trace` when trimming left flies with
different lengths.

Without this, `OverlapsSleep`/`OverlapsWake` are `NaN`, the wake cross-tab
is empty, and both `validate_shuffling` checks skip with a reason.

### 3. Upstream death flags

```matlab
opts.DiedMidRecording = logical(64x1);   % from process_LocomToSleep
```

No parallel dead-fly detector is built here; that decision stays upstream.
Flies marked this way are **not** excluded — their pre-truncation data is
real and `Window` already restricts the analysis to it. `shuffling_stats`
pools them by default; pass `IncludeStatuses = {'OK'}` for the strict
reading.

## Two decisions worth knowing about

**Core vs. outer duration.** An episode's first and last runs are bounded on
their outer side by something that is *not* alternation. A fly that sits on
beam 15 for two hours, shuffles `15,14,15`, then sits for another hour
produces one 3-run episode whose outer span is three hours. `DurationSec`
reports that span (it is what `StartSample`/`EndSample` delimit, and the
overlap check needs it); **`CoreDurationSec` excludes both boundary runs and
is the honest measure**. Statistics and the `ChatteringArtifact` threshold
use core.

**Episodes may share a boundary run.** After closing an episode the scan
resumes *at* its last run, so `[15 16 15 14 15 14]` yields a 3-run and a
4-run episode sharing run 3. That is deliberate — the shared run really does
participate in both patterns — but it means **episode durations must never
be summed** to get total shuffling time. Coverage masks are used throughout
instead; `PerFly.ShuffleSeconds` is already deduplicated.

## Deferred

Day/night and any ZT split is **not** computed. It needs the phase mask
`segment_sleep_phases` builds internally but does not return, and the
`.light`/`.dark` fields cannot substitute — they are non-contiguous
concatenations, so indexing them by sample would give the wrong phase.
`stats.DayNightIncidence` carries this statement so the gap is visible in
the output rather than silently absent. When added, pair it per fly
(one day rate and one night rate each) rather than pooling — the design is
repeated measures.

## Testing

```matlab
test_shuffling_detection
```

22 tests, no toolboxes required. Test 15 is the one that matters most: a
synthetic fly with a known episode at a known time, checked to recover the
same elapsed time *and* the same sleep bin after a 1200-sample trim. It has
been mutation-checked — dropping the offset in either place makes it fail.

Written against the MATLAB/Octave common subset (no `table`, no `string`, no
`inputParser`), so it runs on older MATLAB. Verified under Octave 8.4;
`PreProcessV2`'s cleaning was checked against `PreProcess`'s exact loop on
randomized inputs.
