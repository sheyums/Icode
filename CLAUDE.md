# Icode — left-truncated bout-duration fitting

MATLAB (R2025a+) maximum-likelihood fitters for *Drosophila* sleep and wake
bout durations, plus a model-comparison pipeline. Everything is
**left-truncated**: bouts below a protocol threshold (`xmin`, e.g. 300 s for
the 5-minute sleep criterion) do not exist in the data, so every likelihood
conditions on `T >= xmin`.

## Layout

Two unrelated kinds of fitter live here, and they fit different physical
quantities. The engine and its `Fit*MLE` wrappers fit **bout durations** --
left-truncated positive times. `shiftlognormal_MLE` and
`GeneralizedHyperbolic_MLE` fit **noise distributions** -- the recorded signal
itself, real-valued and untruncated. A log-likelihood from one side is
meaningless against one from the other, and not for the usual reason: it is
not that the criteria penalise differently or the dominating measure differs,
it is that the two are fitted to DIFFERENT DATA. Never let a noise fitter into
a bout comparison, whatever its AIC.

### Bout-duration suite

| File | Role |
| --- | --- |
| `FitTruncatedDiscreteMLE.m` | **Shared engine.** One truncated likelihood, nine families in a `switch` registry. Adding a family costs ~15 lines: `ParamNames`, `cdf`, `sf`, `logpdf`, `unpack`, `valid`, `starts`; a mixture also needs `report`/`nReport` and a `guard`. |
| `FitHyperexponentialMLE.m` | Mixture of exponentials, K=1..N, parametrized in **observed** weights `q`. |
| `FitExponentiatedWeibullMLE.m` | Exponentiated Weibull; `FixAlpha=true` nests the plain Weibull. |
| `FitHyperErlangMLE.m` | Mixture of Erlangs, shapes `[1..1 m]` swept over m. The exponential-native route to a **non-monotone hazard**. The sweep is warm-started (`WarmStart`, `SweepStarts`) via the engine's `StartZ`; cold starts are reduced, never removed, so a warm start can only add a candidate optimum. Refuses rather than return a fit whose guard rejected every shape. |
| `FitWeibullMixtureMLE.m` | Two-component Weibull mixture; the other non-monotone-hazard family. |
| `FitGammaMLE` … `FitBetaMLE` (8 files) | Thin wrappers over the engine. |
| `CompareBoutModels.m` | Fits the whole library, ranks by AICc/BIC, G-tests the winner, plots it. `pearson3` and `beta` are **opt-in**, not part of a default run. |
| `HyperexponentialLRT.m` | Parametric bootstrap LRT for mixture order. |
| `BoutHazard.m` | Life-table hazard with Wilson bands. Does NOT rank models — it decides whether a whole FAMILY can work, by asking whether the hazard rises. |
| `test_*.m` (6 files) | 129 tests. Run each by name from this directory **in MATLAB**. See Testing for the Octave caveat. |

**The bout-duration library is CLOSED.** The families `CompareBoutModels`
fits are the only ones to be used for bout distributions -- the user's
decision, 2026-09-14 local. By default: `hyperexponential` (K=1..N),
`exp_weibull`, `weibull`, `gamma`, `powerlaw_cutoff`, `powerlaw`, `erlang`,
`chisquared`, `weibull_mix`, `hyper_erlang`; `beta` and `pearson3` opt-in by
name. The engine's registry is cheap to extend (see its row above) but that is
a statement about the CODE, not an invitation: do not add a family, and do not
port one in from elsewhere, unless asked for it directly.

The analysis folder holds about 37 further `*_MLE` files from roughly ten years
of separate work, and **most of them DO fit durations** -- `ExpoMLE`,
`GammaMLE`, `LogNormMLE`, the `*ExpoMLE` mixtures, `PlMLE` and others. (An
earlier version of this file said they fit noise or other quantities. That was
wrong; only `shiftlognormal_MLE` and `GeneralizedHyperbolic_MLE`, in the
section below, are noise fitters.) They are not candidates for three reasons
that have nothing to do with what quantity they fit: the user decided the
library is closed, they are **continuous** where `CompareBoutModels` is
discrete by default, and most condition on **`[xmin, xmax]` with `xmax`
defaulting to `max(data)`** -- double truncation whose normalizer depends on
the largest observation, which is a different likelihood, not a different
parametrization of ours.

If a dataset cannot be described by this library, that is a FINDING to report
-- as the per0 wake bouts were, where every monotone family failed and the
hazard said why -- not a cue to widen the library.

### Noise-distribution fitters — NOT bout fitters

These fit the noise in the recorded signal, not durations. Nothing in
`CompareBoutModels`, `BoutHazard` or the engine applies to them: no `xmin`, no
`SamplingInterval`, no truncation, no hazard argument, and no shared AIC table.

| File | Role |
| --- | --- |
| `shiftlognormal_MLE.m` | Shifted-lognormal noise fitter. Positive support with a free shift, **untruncated**. |
| `GeneralizedHyperbolic_MLE.m` | GH / NIG noise fitter, five physical parameters (`mu`, `lambda`, `alpha`, `beta`, `delta`); `AUTO` chooses GH vs NIG by LRT and BIC. **Real-valued support and untruncated** -- the whole real line, so it does not even share a support with a duration law. **Open, unfixed:** on a GH fit to data from the `delta -> 0` variance-gamma limit, `profile_ci.delta` returned `[2.417e-305, 0.9076]` -- an underflowed finite number where the header promises 0 or Inf for an open bound, apparently the root search walking out to `log(delta) ~ -700`. Reported by another session; the user has seen it and has not asked for a fix. Do not trust a `profile_ci` bound that is subnormal. |

### Signal-level tools — fit no distribution at all

| File | Role |
| --- | --- |
| `chi2p.m` | Sokolove-Bushell chi-square periodogram against a block-permutation null. Makes no assumption about waveform shape, so a sharply peaked circadian profile scores on equal footing with a sinusoid. Answers *what period*, not *what distribution*. |
| `jsd_kde.m` | Jensen-Shannon distance between two samples by KDE, with bootstrap CIs and a noise-floor correction. Compares two empirical distributions; fits neither. Calls `ksdensity` (Statistics Toolbox). |
| `stressTest_jsd_kde.m`, `stressTest_chi2p.m` | Stress suites for those two: reported 36 and 44 tests (measured by another session in MATLAB R2026a at 37cf612; not re-measured here). **Not independent validation** -- per the git history in the analysis repo, the suites AND the code they test were both Claude-authored (every commit touching `chi2p`, the suites, and the recent `jsd_kde` / `GeneralizedHyperbolic_MLE` edits carries a Claude co-author line). A passing suite here means self-consistency, not a second opinion. Named `stressTest_*`, so anything globbing `test_*` misses them. |

## Analysing a new dataset

Works for any left-truncated durations or sizes, not only sleep and wake
bouts. **The order matters**: the hazard decides which families are
ADMISSIBLE, the likelihood decides which admissible one is BEST, and the
guards decide whether the winner's parameters mean anything. Run it the other
way round and you get a confident fit from a family the data had already
excluded.

**0. Establish `xmin` and `dt` first. They are facts about the protocol, not
choices.** `xmin` is the threshold below which durations do not exist -- 300 s
for a 5-minute sleep criterion, or just the scoring resolution if nothing else
censors. `dt` (`SamplingInterval`) is the recording grid. Everything conditions
on `T >= xmin` and bins on `dt`; a wrong `dt` silently destroys resolution, and
the grid-mismatch guard will name the spacing it detects instead.

**1. Look at the hazard BEFORE fitting anything.**

```matlab
H = BoutHazard(x, xmin, SamplingInterval=dt, NumBins=16, Plot=true);
```

A hazard that falls and then RISES excludes -- by construction, not by
evidence -- every hyperexponential at every order (a sum of decreasing
exponentials is decreasing), every other monotone family here, and the
exponentiated Weibull (one turning point at most). It also rules out
between-individual heterogeneity as the explanation, since mixing shifts weight
toward longer-lived components and can only make a hazard fall FASTER. No
information criterion will tell you any of this: they rank the candidates you
thought of. For a test rather than a picture, pass `NullSurvivalHandle` (a
fitted hyperexponential's `SurvivalHandle`) and `NullReplicates`, fix the
`NumBins` set before looking, and check `H.NullDropped` is 0.

**2. Fit the library.**

```matlab
R = CompareBoutModels(x, xmin, SamplingInterval=dt, MaxComponents=5);
```

Read it in this order: **`R.Skipped`** (a model that never competed is not a
model that lost), then `R.Table` (ranking, `Degenerate`, `Reason`), then
`R.Selected` and `R.SelectionPath` -- the walk-down stops at the first row that
passes the G test, so rows below the winner are UNTESTED, not accepted.

**3. Check the winner before believing it.**

| field | what it tells you |
| --- | --- |
| `Diagnostics.GuardOK` | false means parameters are unidentified, whatever the likelihood says |
| `Diagnostics.MixtureMinCount` | observations in the smallest branch; a handful means that branch is decoration |
| `Diagnostics.TailFraction` | `S(xmin)`, how much of the untruncated law survived truncation |
| `q` vs `Theta`'s `w` | report `q`. `w` is extrapolation below `xmin` |
| branch MEANS, not stage counts | `m/lambda` is identified; the integer `m` usually is not -- read `ShapeSweep` and see how flat it is |
| `k` vs `numel(Params)` | criteria use `k`; mixtures print one number more than they charge |

**4. When two orders are close, test rather than rank.** AICc and BIC can
disagree by design (they penalise differently), and neither answers "is this
component real". `HyperexponentialLRT` does, with a simulated null -- the
mixture-order LRT is non-regular, so no chi-square applies. `CompareBoutModels`
triggers it automatically when either criterion separates neighbouring orders
by less than `OrderLRTThreshold`.

**5. Prefer a comparison at MATCHED k.** When two candidates have the same
free-parameter count, the information-criterion penalties cancel and the
difference is pure likelihood -- immune to every argument about how `k` is
counted. On the per0 wake bouts, `hyper_erlang` against `hyperexp K=3` (both
k=5) differed by 18.7 nats, and `dAICc = 37.41` is exactly twice that.

**6. State the caveats that the numbers cannot.** If a family was added
BECAUSE the earlier library failed on these data, its goodness-of-fit p-value
is optimistic -- the family was chosen having seen the data. Say so, or re-test
on held-out individuals. Structural arguments (what a family's hazard CAN do)
do not carry that problem and are the stronger claim.

## Conventions that are not optional

**`SamplingInterval` is required everywhere**, in the same units as `xmin`.
It sets the discrete bin width and `n_min = max(1, round(xmin/dt))`. A wrong
value silently destroys resolution; the grid-mismatch guard warns and names
the spacing it detects in the data.

**`DistributionType` defaults to `"discrete"`.** Discrete log-likelihoods are
unit-free; continuous ones shift by `n*log(dt)`. **Never compare AIC/AICc/BIC
across modes** — different dominating measure, not different evidence.

**Report `q`, not `w`.** `q_j` is the observed-population weight (the fraction
of retained bouts from component *j*); `w_j` is the untruncated weight, related
by a factor `exp(+lambda_j * xmin)` that can reach 240x. When a time constant
sits below `xmin`, `w` is extrapolation into a region the protocol excluded.
Check `WeightsUntruncatedReliable`.

**`k` is not `numel(Params)`.** `Params` reports every parameter; `k` counts
the free ones. Hyperexponential at K=3 reports 6 numbers but has k=5, because
`sum(q)=1`. Information criteria use `k`. `nReported` records the difference.

**`q` is the OBSERVED weight in EVERY mixture** — hyperexponential, `hyper_erlang`
and `weibull_mix` alike. It is component *j*'s share of the bouts that were
actually retained, `q_j ~ w_j S_j(xmin)`. The fitted weights inside
`hyper_erlang` and `weibull_mix` are untruncated `w`; `M.report` converts them
at output. Never print a `w` under the name `q`: one fit had `w = 0.849` on a
component holding `2.5e-18` of the data.

**`Theta` is not `Params`.** `H.Theta` is the model's natural parameter vector,
what `M.sf`/`M.cdf`/`M.guard` are written in and the only thing safe to feed
back into the model. `H.Params` is what gets *reported*, and for a mixture it is
a transform of `Theta` with its own delta-method SEs. `nReport` exceeds `nPar`
there, because `sum(q)=1` leaves one weight determined but still worth printing.

**`pearson3` is opt-in: name it in `Models=` or it is not fitted.** It walked a
ridge to the optimizer's evaluation budget on every start -- 516 s for THREE
starts on 1500 bouts, so ~70 min for the engine's 26 -- and its own guard then
rejected the result, so the default library was spending most of its compute on
a row that could not be selected. `Models={'pearson3'}` still fits it, measured
at 35 s for n=120 in Octave, and `Degenerate` came back 0 there: the ridge
needs enough data to become the attractor, so do not read the guard as
unconditional. Fit it deliberately when you want to see the ridge -- on 1500
bouts the REJECTED fit outscored every legitimate model (`FINDINGS.md`).

**A model absent from `R.Table` did not necessarily lose.** `R.Skipped` lists
candidates that never entered the comparison, with the reason — usually a
mixture whose guard rejected every configuration, which is a statement about
the data rather than a malfunction. Absent from both is a bug; absent from the
table but present in `R.Skipped` is a finding. One case is in BOTH on purpose:
a `hyper_erlang` refused at `HyperErlangComponents=J` steps down to J-1..2,
so the refusal stays in `R.Skipped` while the first identified order enters
the table as `hyper_erlang J=j` — named for the order it was fitted at, never
the bare family name. `HyperErlangStepDown=false` turns it off.

**`S(xmin)`** is the *untruncated* survival at the threshold — the likelihood's
normalizer, reported as `Diagnostics.TailFraction`. **`SurvivalHandle`** is the
*truncated* survival `S(t)/S(xmin)`, which is 1 at `xmin`. Different objects;
one is a scalar, the other a curve.

**`BoutHazard`: only `RiseNullP` is a test, and it depends on `NumBins`.**
`RiseRatio` is a maximum taken after a minimum, so it is >= 1 for any curve;
`RiseCI` bootstraps the data, not a null, and is DESCRIPTIVE only (its lower
limit exceeded 1 on 9 of 50 no-rise datasets). A test needs
`NullSurvivalHandle` (a fitted monotone model, e.g. a hyperexponential's
`SurvivalHandle`) and `NullReplicates`; without them `NonMonotone` falls back
to `RiseDisjoint`, which is pointwise. `RiseNullP` moves with the bin count --
0.010 at 16 bins and 0.068 at 20 on the same wake bouts -- so **fix the set of
`NumBins` before looking, and report every one** (FINDINGS judgement call 8).

## Traps that have already bitten, with tests guarding them

- **A mixture guard must test the OBSERVED share, not the mixing weight.**
  `min(w)*n` misses the failure it exists to catch: under left truncation a
  component whose mass sits below `xmin` carries a large `w` and holds none of
  the data. A 3-component fit at `xmin=100` had `w = [0.055, 0.809, 0.136]`
  with rates `[0.019, 0.566, 0.0074]` — the LARGEST weight had mean 1.8 s, so
  `S(100)=3e-25` and it held zero of 1500 observations. `min(w)*n` read 82.3.
  Both `heGuard` and `weibullMixGuard` now test `w_j S_j(xmin)`, normalized.
  The phantom component had bought a stage count the data never supported.
- **The same truncated law has several parametrizations, and only some look
  degenerate.** The fit above and one with `w = [0.288, 1.3e-10, 0.712]` have
  identical log-likelihoods to ten digits — they ARE the same distribution
  above `xmin`. Judge identifiability on what the data see, never on where the
  optimizer happened to stop.
- **Never pass a name-value pair twice.** Octave's generated twins parse
  `varargin` in a loop and take the last; MATLAB's `arguments` block rejects a
  duplicate name outright. A one-platform failure the twins cannot catch,
  because the twins are what replace the `arguments` block.
- **Integer-shape Erlang needs no `gammainc`.** `S = exp(-x) sum_{j<m} x^j/j!`
  is exact, ~11x faster for the survival and ~6x for the CDF, and more accurate
  in the deep tail than Octave's `gammainc` (which drifts to 1.5e-3 relative at
  `x=0.1, m=8`).
- **Benchmark the function the hot path CALLS, not its sibling.** That 11x was
  measured on `erlangSFint`. The discrete likelihood calls the *CDF* for every
  bin edge, and `erlangCDFint` went unbenchmarked: it summed the upper tail
  directly whenever `x < m`, which is most of the fitted range for a slow
  branch, and was **0.7x** — 40% SLOWER than the `gammainc` it replaced. A
  full-library call that had cost 379 s was still running after 17 minutes,
  while the commit message advertised a speedup.
- **Order numerical branches by what threatens the answer, not by the obvious
  split.** `F = 1 - S` cancels with relative error about `eps/F`: harmless at
  `F = 0.1`, fatal at `F = 1e-300`. Branching on `x < m` put the slow path in
  the common case; branching on `F < 1e-6` puts it only where cancellation
  actually bites — and those are the small-`x` elements where the tail
  converges in a few terms, so the accurate branch is also the rare and cheap
  one.
- **The same seed is NOT the same data across implementations.** `rng(25)`
  then `rand` gives different draws in MATLAB and Octave, and a generator that
  rejection-samples diverges on the first rejection. So a log-likelihood
  measured in one and compared against the other is comparing two datasets, not
  two fitters. This produced a "77-nat" gap that was really 13 once both sides
  were measured on the same sample. Cross-platform claims need the same numbers
  computed twice, or the data written to a file and read by both.
- **Octave is not a proxy for MATLAB.** `*_oct.m` twins are local scaffolding
  (gitignored, regenerable) used for development. Octave has accepted three
  things MATLAB rejects: growing a struct array from `struct([])`,
  `struct2table` on a 1x1 struct with a `''` field, and it hid a false-positive
  guard because output was filtered to `[FAIL]` lines. **Anything touching
  struct arrays, tables or graphics needs a MATLAB run.**
- **A new MATLAB figure's theme is not predictable — pin it.** In R2026a
  batch mode, `figure('Name',..)` came out light while a `'Visible','off'`
  figure a moment later came out dark (axes `Color` 0.07, `XColor` 0.85), and
  another session got the opposite. `CompareBoutModels`' figure for the per0
  wake bouts rendered dark and its black data curve was invisible.
  `figure('Color','w')` alone does NOT fix it: the margin goes white while the
  axes stay dark. Create the figure, then `try theme(fh,'light'); catch, end`
  BEFORE any axes, so children inherit (`theme` is absent in Octave and before
  R2025a). `plotFit` does this and test 23 asserts a white figure and axes.
- **A statistical criterion is not validated by tests on single datasets —
  calibrate it on many under a null.** `BoutHazard`'s first rise criterion
  passed all 10 of its tests, including known-hazard laws, because they asked
  "does it detect a rise that exists, and decline one that doesn't" on one
  dataset each. Its error rate was never measured. Run on 50 datasets from a
  strictly decreasing hazard, it claimed a rise in 16-20% of them. Before
  quoting any p-value or accept/reject rule, simulate from the null it claims
  and count how often it fires.
- **A null replicate must pass through every data-dependent step the observed
  statistic did.** `RiseNullP` first binned each null dataset on the OBSERVED
  data's edges and at-risk/saturation mask, while the observed value had been
  computed after choosing those from itself. Not exchangeable, and it pushed p
  from 0.019 to 0.048 (16 bins) and from 0.057 to 0.173 (20 bins). A sampler one
  grid step early, found in the same code, moved p by at most 0.012 there --
  within Monte Carlo error -- on the real bouts with a fitted null. That is not
  a general reprieve for discretization (next entry): measure which defect
  matters IN THE SETTING AT HAND before fixing the one that looks worse. Fixed
  in 2371913: each replicate runs the whole estimator on its own edges, and the
  sampler uses `pmf(k) = S(g_k - dt) - S(g_k)`.
- **A calibration harness must discretize its simulated data exactly as the
  estimator's null sampler does.** Calibrating `RiseNullP` (pre-2371913) on
  synthetic data generated with `round()` -- nearest bin -- while the null
  sampler followed the engine's convention (a bout of k bins lies in
  `((k-1)dt, k dt]`, i.e. `ceil`) gave min p 0.080, median 0.610, 0/40 below
  0.05: a clean, confident "the test is conservative". Changing only the
  harness to `ceil()` gave min 0.025, median 0.542, 3/40 below 0.05 --
  approximately uniform. The conservatism was the HARNESS. The two traps above
  are about the estimator; this one is about the thing it is tested with,
  which nothing else here guards. The same half-bin question barely moved p on
  the real bouts, so which discretization defect dominates depends on the
  setting (a plausible, unmeasured reason: log-spaced edges follow `max(t)`,
  which moves little across synthetic draws and a lot between a real sample
  and a null draw). Simulate through the engine's own survival convention --
  `sampleFromSurvival`-style inversion of `SurvivalHandle` -- not `round()`.
- **Never floor a probability without its normalizer.** Flooring a bin
  probability at `realmin` while `S(xmin)` sat in the subnormals made `p/S`
  reach 4.5e15 — a conditional probability above 1, worth +36 log-likelihood
  per observation. A Weibull scored `logL = +39972` and won a comparison.
  Engine test 23 asserts `logL < 0`.
- **A small `S(xmin)` is not a broken fit.** For a family anchored at 0 it only
  means the untruncated interpretation is vacuous. A gamma at `shape -> 0` has
  `S(xmin) ~ 3e-12` and is sound: `Gamma(a,z) -> E1(z)`, so the truncated law
  tends to the power-law-with-cutoff limit. Gating on it excluded good fits.
  Only a **free location** produces a real ridge — see `pearson3Guard`.
- **Don't write tests that assert where the optimizer lands.** These
  likelihoods are non-concave; MATLAB and Octave find different legitimate
  optima. Assert the decision *rule* instead (engine test 24).
- **Some rows are the same model twice.** `powerlaw_cutoff` IS `gamma`
  (`alpha = 1-shape`, bit-identical logL); `erlang` at shape 1 IS `hyperexp K=1`
  IS the exponential. `R.Nesting` lists them. Equal log-likelihoods there are an
  implementation check, not a coincidence — disagreement beyond ~1e-9 is a bug.
- **`chisquared` is not unit-invariant** (its scale is fixed at 2), so its fit
  depends on whether you scored in seconds or minutes. Exclude it rather than
  interpret its defeat.

## Testing

```matlab
test_FitTruncatedDiscreteMLE      % 28
test_FitHyperexponentialMLE       % 26
test_FitExponentiatedWeibullMLE   % 23
test_CompareBoutModels            % 29
test_HyperexponentialLRT          % 13
test_BoutHazard                   % 10
```

**In MATLAB, run them by name from this directory.** Check `which
test_CompareBoutModels -all` first: a stale copy elsewhere on the path wins
over the repo unless you are `cd`'d into it, and that has already caused a
"pass" that ran 24 of 27 tests.

**Anything that moves `pwd` away from this repo defeats that check**, and
`run('some\path\script.m')` does exactly that: MATLAB changes into the
script's folder for the duration. The mechanism is NOT that the script sits
beside stale copies — the scratch folder may contain no fitters at all. It is
that leaving the repo removes the current-folder precedence that was making the
repo win, so resolution falls through to the MATLAB **path**, and whatever
other copy sits on the path takes over. `which -all` run from the repo
therefore reports the repo copy winning while `run()` quietly uses another one.
This has already produced a discarded 35-minute timing measured against a
pre-fix engine.

Either `addpath` the scratch folder and call the script by NAME, or open with

    assert(isequal(fileparts(which('CompareBoutModels')), pwd))

so a run against the wrong copy fails loudly instead of returning numbers.

**Do not "clean up" the other copies.** The folder that shadows this one is the
user's main analysis repo, not a pile of stale files: its
`FitHyperexponentialMLE.m` is tracked there and called by scripts in it, and
its copies of the newer fitters are UNTRACKED, so deleting them cannot be
undone. Fix the resolution, not the files — take that folder off the path while
working here, or refresh its copy of whatever you need.

**That folder now carries `Icode_sync_stamp.txt`**, naming the commit its
copies were synced from (`1397554`, 2026-09-14) with a SHA-256 prefix per file
and the backup location of what was replaced. A session working from there
should compare that file's Commit line against `git log -1` in Icode BEFORE
trusting the copies. And note the asymmetry: `addpath('...Icode','-begin')`
does NOT protect you, because MATLAB's current folder beats the path — the
user's analysis scripts live in that folder, so its copies win whenever it is
`pwd`. The only durable options are to re-sync after any `.m` change, or to
keep one copy of each function.

**In Octave they do NOT run by name.** Three of the suites call the bare
`Fit*MLE` names, which resolve to the MATLAB originals and die on the
`arguments` block. They need a generated test twin — see `verify_all.sh` in the
scratchpad. `test_CompareBoutModels` dispatches to `_oct` names itself and does
run directly.

No toolbox required — local Marsaglia-Tsang gamma sampler, and chi-square tails
via `gammainc(...,'upper')` rather than `chi2cdf`. Engine test 20 cross-checks
against `gamcdf`/`wblcdf` and self-skips when absent, so it only runs in MATLAB.
