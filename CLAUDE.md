# Icode — left-truncated bout-duration fitting

MATLAB (R2025a+) maximum-likelihood fitters for *Drosophila* sleep and wake
bout durations, plus a model-comparison pipeline. Everything is
**left-truncated**: bouts below a protocol threshold (`xmin`, e.g. 300 s for
the 5-minute sleep criterion) do not exist in the data, so every likelihood
conditions on `T >= xmin`.

## Layout

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
| `test_*.m` (5 files) | 119 tests. Run each by name from this directory **in MATLAB**. See Testing for the Octave caveat. |
| `shiftlognormal_MLE.m` | Pre-existing noise fitter. **Untruncated** — do not put it in an AIC table with the others. |

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

**In Octave they do NOT run by name.** Three of the suites call the bare
`Fit*MLE` names, which resolve to the MATLAB originals and die on the
`arguments` block. They need a generated test twin — see `verify_all.sh` in the
scratchpad. `test_CompareBoutModels` dispatches to `_oct` names itself and does
run directly.

No toolbox required — local Marsaglia-Tsang gamma sampler, and chi-square tails
via `gammainc(...,'upper')` rather than `chi2cdf`. Engine test 20 cross-checks
against `gamcdf`/`wblcdf` and self-skips when absent, so it only runs in MATLAB.
