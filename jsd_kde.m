function [d, d_CI, d_SE, d_boot, noise_P, noise_Q, d_corr, d_corr_CI, d_corr_SE, p_perm, d_perm] = jsd_kde(P, Q, Ngrid, min_cutoff, do_log_transform, n_boot, ci_alpha, rng_seed, n_perm, Verbose)
% JSD_KDE  Jensen-Shannon distance between two univariate distributions
%          estimated via kernel density estimation (KDE).
%
% SYNTAX
%   d = jsd_kde(P, Q)
%   [d, d_CI, d_SE, d_boot, noise_P, noise_Q, d_corr, d_corr_CI, d_corr_SE] = ...
%       jsd_kde(P, Q, Ngrid, min_cutoff, do_log_transform, n_boot, ci_alpha, rng_seed)
%   [..., p_perm, d_perm] = jsd_kde(..., rng_seed, n_perm)
%   ... = jsd_kde(..., n_perm, Verbose)
%
% DESCRIPTION
%   Estimates the Jensen-Shannon distance (JSD; the square root of the
%   Jensen-Shannon divergence) between the empirical distributions of two
%   sample vectors P and Q.  Continuous density estimates are obtained via KDE
%   with an automatically selected bandwidth (ksdensity's normal-reference
%   plug-in rule), evaluated on a fixed grid.  Numerical integration uses the
%   trapezoidal rule.
%
%   For right-skewed data a log-transform is optionally applied before KDE.
%   The divergence is then evaluated in the transformed space directly; see
%   SCALE INVARIANCE below for why this needs no Jacobian correction.
%
%   When n_boot > 0, a nonparametric bootstrap yields confidence intervals and
%   a noise-corrected distance, d_corr.  Finite-sample KDE bias is estimated by
%   computing JSD between pairs of bootstrap resamples drawn from the *same*
%   distribution (P vs P', Q vs Q').  The mean bias is taken as an ensemble
%   average over all n_boot replicates and subtracted as a constant scalar in
%   divergence (squared-distance) space before the square root.
%
%   When n_perm > 0, a label-permutation test reports p_perm, the probability
%   of a distance at least this large when P and Q are relabellings of one
%   pooled sample.  This is the only output of this function that constitutes a
%   hypothesis test; the bootstrap quantifies and de-biases d but cannot say
%   whether it exceeds chance.  See PERMUTATION NULL.
%
% INPUTS
%   P, Q              Numeric column (or row) vectors of observations.
%                     All values must be >= min_cutoff.  Neither may be empty.
%   Ngrid             Number of evaluation grid points (default: 512).
%   min_cutoff        Hard lower bound of the support.  Pass [], NaN, or omit
%                     it to take the default of min([P; Q]).  Set to a known
%                     physical minimum (e.g. 0 for positive-definite
%                     quantities) to enable bounded-support
%                     KDE with boundary correction.  NOTE: at its default the
%                     cutoff coincides with the observed minimum, which selects
%                     the unbounded estimator; bounded support is used only
%                     when min_cutoff lies strictly below the data.
%   do_log_transform  Logical scalar; if true, log-transforms the data before
%                     KDE.  If empty or omitted, auto-detected from sample
%                     skewness (default).  Auto-detection fires on default
%                     calls, so a skewed sample takes the log path unless the
%                     transform is explicitly disabled with false.
%   n_boot            Number of bootstrap replicates (default: 0, no bootstrap).
%                     The replicate loop is serial; see ALGORITHM NOTES for
%                     why it is not parallelised.
%                     For publication-quality CI and noise correction,
%                     1000-2000 is a reasonable choice.
%   ci_alpha          Nominal coverage level, alpha in (0,1)
%                     (default: 0.05 => 95% CI).
%   rng_seed          Nonnegative integer seed for the bootstrap, or [] to draw
%                     from the current stream (default: []).  See
%                     REPRODUCIBILITY below.
%   n_perm            Number of label-permutation replicates (default: 0, no
%                     permutation test).  Independent of n_boot: the bootstrap
%                     quantifies uncertainty in d, the permutation test asks
%                     whether d exceeds chance.  See PERMUTATION NULL below.
%                     1000-5000 is reasonable; the attainable p-value is
%                     bounded below by 1/(n_perm + 1).
%   Verbose           Logical scalar (default: true).  When false, suppresses
%                     the jsd_kde:logTransform announcement only.  Intended for
%                     callers that invoke this function in a loop, where that
%                     announcement fires on every call and drowns everything
%                     else.  It does NOT silence any other warning: degenerate
%                     input, a skipped permutation test and degenerate
%                     replicate counts all still report, because those bear on
%                     whether the returned numbers mean what they claim.
%
% OUTPUTS
%   d                 Jensen-Shannon distance, scalar in [0, 1] (bits basis).
%                     Deterministic given the data; the bootstrap does not
%                     enter it.
%   d_CI              1x2 bootstrap CI on d; [] if n_boot = 0.
%   d_SE              Bootstrap standard error of d; [] if n_boot = 0.
%   d_boot            n_boot x 1 bootstrap distribution of d; [] if n_boot = 0.
%   noise_P           Mean bootstrap JSD distance of P vs P' (noise floor).
%   noise_Q           Mean bootstrap JSD distance of Q vs Q' (noise floor).
%                     Both are descriptive distance-scale summaries only; the
%                     actual correction is performed in divergence space.  Note
%                     these return 0, not [], when n_boot = 0.
%   d_corr            Noise-corrected JSD distance: ensemble-mean KDE bias
%                     subtracted in divergence space, then sqrt; clamped at 0.
%                     Returns 0 when n_boot = 0, since no noise estimate exists.
%   d_corr_CI         Bootstrap CI on d_corr; [] if n_boot = 0.  NOT a
%                     hypothesis test: d_corr is clamped at 0, so when the true
%                     distance is near zero the bootstrap distribution piles up
%                     against that boundary and the interval loses its nominal
%                     coverage.  "d_corr_CI excludes 0" is therefore not
%                     evidence that P and Q differ; use p_perm for that.
%   d_corr_SE         Bootstrap standard error of d_corr; [] if n_boot = 0.
%                     Same boundary caveat applies -- an SE is a poor summary
%                     of a distribution truncated at zero.
%   p_perm            One-sided permutation p-value for H0: P and Q are drawn
%                     from the same distribution.  NaN if n_perm = 0, and NaN
%                     on the degenerate and zero-range early-return paths,
%                     where no test was performed.
%   d_perm            n_perm x 1 null distribution of the JSD distance under
%                     label permutation; [] if n_perm = 0.
%
% SCALE INVARIANCE
%   JS divergence is exactly invariant under a smooth invertible
%   reparametrisation y = g(x): the Jacobian cancels inside the log because the
%   mixture transforms identically to its components,
%       p_Y/m_Y = (p_X/g')/(m_X/g') = p_X/m_X,
%   and the remaining g' cancels against dy = g' dx.  The estimator therefore
%   evaluates the integrals directly in the transformed space on a uniform grid
%   rather than back-transforming to the original scale.  This is the same
%   estimator either way, but the uniform grid and the absence of a division by
%   exp(y) make the trapezoidal quadrature markedly more accurate.
%
%   A corollary worth knowing: the log-transform offset (`shift`) cannot bias
%   the estimand, only the smoothing, since it is part of a monotone map.
%
% PERMUTATION NULL
%   The bootstrap and the permutation test answer different questions, and
%   neither substitutes for the other.
%
%     noise_P/noise_Q/d_corr  estimate and remove finite-sample KDE BIAS.  Each
%       noise term resamples ONE group against itself, so the reference
%       distribution is that group's own shape.  This is a bias correction; it
%       yields no p-value, and the bootstrap CI on d_corr is a CI on a
%       boundary-clamped statistic, so its lower limit is not a hypothesis test.
%
%     p_perm/d_perm  estimate the NULL DISTRIBUTION of d.  The transformed
%       samples are pooled, then re-split WITHOUT replacement into pseudo-groups
%       of exactly nP and nQ.  Any d this produces is attributable to sampling
%       alone, because the labels carry no information by construction.  The
%       reference distribution is the POOLED shape, which is what H0 asserts
%       both groups share.
%
%   Held fixed across replicates, and this is load-bearing: the log-transform
%   decision, the shift, min_cutoff, min_is_observed and the evaluation grid.
%   Permutation operates on P_trans/Q_trans, i.e. AFTER the transform, which is
%   legitimate because a monotone map applied identically to both samples
%   commutes with relabelling.  Re-deciding the transform per replicate would
%   be a bug: do_log_transform auto-detects from per-sample skewness, so a null
%   built by recursive calls into this function would silently average
%   transformed and untransformed replicates.  Callers writing their own
%   permutation wrapper around jsd_kde must pass do_log_transform and
%   min_cutoff explicitly for exactly this reason.
%
%   Bandwidth is re-estimated per replicate, as it is in the bootstrap loop.
%   The observed d uses bandwidths selected from samples of size nP and nQ, so
%   a faithful null must pay the same bandwidth-selection cost at the same
%   sizes.  Holding bandwidth fixed would understate the null precisely when
%   one group is small, which is the case the test is most needed for.
%
%   What this null does and does not control for.  It controls for sample size,
%   for size ASYMMETRY between the groups, and for the shape of the pooled
%   distribution.  It does NOT need to control for a common rescaling of both
%   groups: JS divergence is exactly invariant to that (see SCALE INVARIANCE),
%   so a uniformly narrower or wider pair of samples has the same estimand.
%   What does move d is the ratio of between-group separation to within-group
%   spread in the TRANSFORMED coordinate, and no null can disentangle that from
%   a genuine effect -- it is the effect.  Reporting the between-group location
%   shift and the within-group spread alongside d is the way to tell a reader
%   which of the two is driving a cross-condition comparison.
%
%   p_perm uses the add-one estimator (1 + #{d_perm >= d}) / (1 + n_perm),
%   which is never zero and stays valid under a finite replicate count.
%
% ALGORITHM NOTES
%   Units.  JSD is computed in bits (log base 2), giving JS divergence in [0,1]
%   and hence d in [0,1].
%
%   Grid extension.  The evaluation grid extends 3 KDE bandwidths beyond the
%   observed data extremes to prevent truncation of KDE tail mass at the
%   integration boundary.  The bandwidth used for this extension is taken on
%   the grid's own scale.  That distinction matters: bounded-support ksdensity
%   works internally on log(x - min_cutoff) and returns a bandwidth on THAT
%   scale, which is meaningless added to raw data values, so a separate
%   raw-scale bandwidth is estimated for sizing in that branch.
%
%   Boundary handling.  Under bounded support the KDE is singular at exactly
%   x = min_cutoff (internal log(0) -> -Inf, then 0/0).  The first grid node is
%   therefore placed just inside the support.  The estimated density tends to 0
%   at the boundary, so no mass is lost.
%
%   Bandwidth.  bw_P and bw_Q are extracted once from the observed samples and
%   used for the point estimate d, for sizing the grid extension, and as the
%   fallback scale on paths where a resample carries no scale information.
%   Inside the bootstrap and permutation loops, however, bandwidth is
%   RE-ESTIMATED from each replicate's own data.
%
%   This was not always so, and the change is worth recording.  Bandwidth used
%   to be cached across replicates, which saved O(6 * n_boot) ksdensity calls
%   but meant no replicate ever paid the cost of choosing its own bandwidth.
%   The bootstrap therefore did not propagate bandwidth-selection variability,
%   leaving d_CI and d_SE narrower than a fully nonparametric bootstrap, and --
%   more seriously -- leaving the noise floor too low, because bandwidth
%   selection is the dominant error source exactly when a sample is small.
%   d_corr consequently under-corrected in the regime it was most needed for:
%   two samples drawn from the SAME lognormal at nP/nQ = 400/10 returned
%   d_corr = 0.198 against a true value of 0, about a third of the raw
%   inflation removed.  Re-estimation closes most of that gap at a real
%   runtime cost, which is the trade this function now takes.  Set n_boot
%   modestly (a few hundred) when scanning many pairs.
%
%   A corollary of re-estimation: a resample or pseudo-group can come out
%   constant even though the observed sample is not, and KDE is undefined
%   there.  Such replicates are scored by the same closed-form singular and
%   identical limits the degenerate-input path uses, and counted in a warning,
%   rather than letting ksdensity improvise a zero bandwidth.
%
%   The bootstrap loop is SERIAL, and deliberately so.  Its six randi calls
%   draw from the global stream in a fixed order, and that draw sequence is the
%   reproducibility contract: it is what lets a saved d_CI be regenerated years
%   later from the same data and rng_seed.
%
%   PARALLELISING THIS COSTS THAT CONTRACT.  A parfor version was written,
%   tested and reverted.  It worked and was correct -- one independent threefry
%   substream per replicate, created on the client and indexed by b, giving
%   results identical across pool sizes -- but substreams necessarily draw
%   DIFFERENT numbers from the serial sequence, so every previously saved
%   d_CI, d_SE, d_boot, d_corr, d_corr_CI and d_corr_SE became unreproducible
%   at the same seed.  Measured on n_boot = 1000: d was bit-identical (it is
%   deterministic and the bootstrap never touches it) while the CI bounds moved
%   by 0.02 and 1.65 standard deviations of the seed-to-seed spread, about as
%   much as simply re-running at another seed.  Statistically equivalent,
%   numerically different -- and for saved results, numerically different is
%   what matters.  The measured speedup was only about 1.5x on two workers.
%
%   If that trade ever becomes worthwhile (many fly-pairs in a loop, say), the
%   working recipe is:
%
%       streams = RandStream.create('threefry', 'NumStreams', n_boot, ...
%                                   'Seed', rng_seed, 'CellOutput', true);
%       parfor b = 1:n_boot
%           idx = randi(streams{b}, nP, nP, 1);   % etc, per replicate
%       end
%
%   What must NOT be done is a bare parfor over the existing randi calls.
%   Seeding the global stream on the client does not seed the workers, each of
%   which has its own, so rng_seed would silently stop reproducing anything and
%   the onCleanup restore below would become meaningless for the draws that
%   matter.  Measured: two parfor runs under an identical rng(42) produced
%   different draws, while the serial loop reproduced exactly.  See chi2p.m,
%   where a closely related mistake -- one RandStream shared across parfor
%   workers, which duplicated surrogates -- was made in earnest and fixed.
%
%   Noise correction.  avg_noise_div is the ensemble mean of all n_boot noise
%   replicates, computed after the loop rather than per-iteration.  This gives a
%   lower-variance noise floor (SE proportional to 1/sqrt(n_boot)) and decouples
%   noise-estimator variance from the corrected bootstrap distribution.
%   Subtracting a scalar shifts that distribution rigidly, so d_boot_corr
%   reflects only the sampling variance of JSD itself.
%
% DEGENERATE INPUT
%   A sample with zero range carries no scale information, so no bandwidth can
%   be estimated from it and the KDE is undefined.  Rather than depend on
%   ksdensity's internal zero-sigma fallback, such input is detected up front
%   and resolved in closed form:
%     - both samples constant at the same value -> d = 0 (identical measures);
%     - both constant at different values, or exactly one constant -> d = 1,
%       since a point mass is mutually singular with respect to any other
%       distribution supported elsewhere, and JS divergence between mutually
%       singular measures is exactly log2(2) = 1 bit.
%   A jsd_kde:degenerateInput warning is issued so the result is never silent.
%
% REPRODUCIBILITY
%   Pass rng_seed (a nonnegative integer) to make the bootstrap deterministic.
%   The bootstrap draws from the global stream, which is seeded on entry and
%   whose previous state is restored on exit via onCleanup, including on error
%   or Ctrl-C, so calling jsd_kde never perturbs the caller's random stream.
%   Omit rng_seed (default []) to draw from the current stream and leave it
%   advanced.
%
%   Only d is deterministic without a seed.  d_CI, d_SE, d_boot, d_corr,
%   d_corr_CI and d_corr_SE all vary run to run, so results intended to be
%   reproducible from saved output should always be generated with a seed.
%   p_perm and d_perm are likewise stochastic.
%
%   The permutation loop runs strictly AFTER the bootstrap loop and draws from
%   the same stream, so every bootstrap output is bit-for-bit unchanged at a
%   given rng_seed whether or not n_perm is supplied.  The converse does not
%   hold: changing n_boot shifts the stream the permutation loop inherits, so
%   d_perm is reproducible only at a fixed (rng_seed, n_boot) pair.
%
%   BREAKING CHANGE, recorded here because a seed alone no longer pins these
%   numbers.  Bandwidth is now re-estimated per replicate (see the Bandwidth
%   note in ALGORITHM NOTES).  d, the resample draw sequence and the evaluation
%   grid are all unchanged, so the resamples themselves are identical to those
%   earlier versions produced -- but d_CI, d_SE, d_boot, noise_P, noise_Q,
%   d_corr, d_corr_CI and d_corr_SE are scored with different bandwidths and so
%   will NOT reproduce values saved by a version that cached bandwidth across
%   replicates.  Expect the noise floor and hence d_corr to move most, and
%   most of all at small or asymmetric sample sizes.  Only d is untouched.
%
% LIMITATIONS
%   d_corr_CI is conditional on the noise floor being known exactly.  Because
%   avg_noise_div is subtracted as a constant, Var[avg_noise_div] is excluded by
%   construction, so the interval is narrower than a full accounting of the
%   uncertainty in d_corr would give.
%
%   When min_cutoff coincides with the observed minimum the estimator is
%   unbounded, so KDE mass falling below min_cutoff is truncated at the
%   integration limit and redistributed by renormalisation.  P and Q are
%   treated identically, so the effect on d is second order.
%
%   Grid resolution is fixed at Ngrid points spanning [min_cutoff, max + ext].
%   A min_cutoff set far below the bulk of the data spends much of the grid on
%   empty space; raise Ngrid in that situation.
%
% REFERENCES
%   Lin, J. (1991). Divergence measures based on the Shannon entropy.
%     IEEE Trans. Inf. Theory, 37(1), 145-151.
%   Endres, D.M. & Schindelin, J.E. (2003). A new metric for probability
%     distributions. IEEE Trans. Inf. Theory, 49(7), 1858-1860.
%   Osterreicher, F. & Vajda, I. (2003). A new class of metric divergences
%     on probability spaces. Ann. Inst. Stat. Math., 55(3), 639-653.
%
% Designed by Sheyum originally but math-checked, re-written and 
% strengthened by Claude and Gemini. January-August 2026.
%

%% ------------------ Defaults & Setup ------------------
arguments
    P           {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeFinite}
    Q           {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeFinite}
    Ngrid       (1,1) {mustBeInteger, mustBePositive}                       = 512
    min_cutoff  double {mustBeScalarOrEmpty}                                = NaN
    do_log_transform                                                        = []
    n_boot      (1,1) {mustBeInteger, mustBeNonnegative}                    = 0
    ci_alpha    (1,1) {mustBeInRange(ci_alpha, 0, 1, "exclusive")}          = 0.05
    rng_seed                                                                = []
    n_perm      (1,1) {mustBeInteger, mustBeNonnegative}                    = 0
    Verbose     (1,1) logical                                               = true
end

P = P(:); Q = Q(:);
nP = numel(P); nQ = numel(Q);

%% ------------------ Reproducibility ------------------
% Seed the global stream when a seed is supplied, and restore whatever state
% the caller had on the way out.  onCleanup fires on normal return, on error
% and on Ctrl-C, so the caller's stream is never left reseeded by this call.
% cleanup_rng must stay in scope for the lifetime of the function; it is
% intentionally never referenced again.
if ~isempty(rng_seed)
    if ~(isnumeric(rng_seed) && isscalar(rng_seed) && isreal(rng_seed) && ...
            isfinite(rng_seed) && rng_seed >= 0 && rng_seed == floor(rng_seed))
        error('jsd_kde:badSeed', ...
            'rng_seed must be [] or a nonnegative integer scalar.');
    end
    rng_state_in = rng;
    cleanup_rng  = onCleanup(@() rng(rng_state_in));
    rng(rng_seed, 'twister');
end

% Resolve data-dependent default. [] and NaN both mean "use the observed
% minimum": [] is what a caller naturally passes to skip a positional
% argument, and the arguments block used to reject it outright.
if isempty(min_cutoff) || isnan(min_cutoff), min_cutoff = min([P; Q]); end

if any(P < min_cutoff) || any(Q < min_cutoff)
    error('jsd_kde:belowCutoff', 'Values fall below min_cutoff.');
end

observed_min = min([P;Q]);
min_is_observed = abs(min_cutoff - observed_min) < 10*eps(observed_min);

%% ------------------ Degenerate (zero-range) input ------------------
% Handled before anything that needs a scale estimate: skewness of a constant
% vector is 0/0 = NaN, and no bandwidth can be estimated from a sample with no
% spread.  Resolving these cases in closed form keeps the result deterministic
% and independent of ksdensity's internal zero-sigma fallback.  Reaching the
% code below also establishes that both samples have positive range, which is
% what guarantees the std fallback for `spread` is strictly positive.
P_is_const = (max(P) - min(P)) <= 10*eps(max(abs(P)));
Q_is_const = (max(Q) - min(Q)) <= 10*eps(max(abs(Q)));

if P_is_const || Q_is_const
    if P_is_const && Q_is_const
        val_tol = 10*eps(max(abs([P(1), Q(1)])));
        if abs(P(1) - Q(1)) <= val_tol
            d = 0;   % identical point masses
            reason = 'both samples are constant at the same value';
        else
            d = 1;   % distinct point masses are mutually singular
            reason = 'both samples are constant at different values';
        end
    else
        d = 1;       % a point mass is singular w.r.t. any spread-out measure
        if P_is_const, reason = 'P is constant'; else, reason = 'Q is constant'; end
    end

    warning('jsd_kde:degenerateInput', ...
        ['Degenerate input (%s): KDE is undefined, returning the exact ' ...
        'limiting distance d = %g.'], reason, d);

    % Every resample of a constant sample is that same constant, so the
    % bootstrap is deterministic and the noise floor is exactly zero.  Outputs
    % follow the same convention as the main path: the noise-corrected
    % quantities are only populated when a bootstrap was actually requested.
    noise_P = 0; noise_Q = 0;
    if n_boot > 0
        d_boot    = repmat(d, n_boot, 1);
        d_CI      = [d, d];
        d_SE      = 0;
        d_corr    = d;
        d_corr_CI = [d, d];
        d_corr_SE = 0;
    else
        d_boot    = [];
        d_CI      = [];
        d_SE      = [];
        d_corr    = 0;
        d_corr_CI = [];
        d_corr_SE = [];
    end

    % No KDE exists on this path, so there is nothing for a permutation null to
    % be a null OF: d came from the closed-form singular/identical limit, not
    % from an estimate.  Report "not tested" rather than manufacturing a
    % p-value, and say so when one was asked for.
    p_perm = NaN; d_perm = [];
    if n_perm > 0
        warning('jsd_kde:permSkippedDegenerate', ...
            ['Permutation test skipped on degenerate input (%s); ' ...
            'p_perm returned as NaN.'], reason);
    end
    return;
end

%% ------------------ Auto-log decision ------------------
% Auto-log functions when min_cutoff is at its default (observed minimum).
if isempty(do_log_transform)
    skew_thresh = 0.8 + 10/sqrt(nP + nQ);
    do_log_transform = (all(P >= min_cutoff) && all(Q >= min_cutoff)) && ...
        ((skewness(P) > skew_thresh) || (skewness(Q) > skew_thresh));
    if do_log_transform && Verbose
        warning('jsd_kde:logTransform', ...
            'Log-transform applied automatically (skewness P=%.2f, Q=%.2f exceeded threshold %.2f).', ...
            skewness(P), skewness(Q), skew_thresh);
    end
elseif do_log_transform && Verbose
    warning('jsd_kde:logTransform', ...
        'Log-transform applied (user-specified).');
end

%% ------------------ Grid Setup (FIXED for Bootstrapping) ------------------
spread = max(median([P;Q]) - min_cutoff, iqr([P;Q]));
% Both samples are known to have positive range here (see the degenerate-input
% block above), so this fallback is strictly positive.
if spread <= 0, spread = std([P;Q]); end

% --- PATCH (3): conditional log-transform shift ---------------------------
% The shift exists solely to keep log(x - min_cutoff + shift) finite and
% well-scaled for the smallest observation.  Let
%     gap = observed_min - min_cutoff   (>= 0)
% be the room already available.  If gap is at least shift_floor no offset is
% needed at all and shift is exactly zero, so the transform is a clean
% log(x - min_cutoff).  Otherwise the offset tops the gap up to shift_floor,
% bounding the smallest transformed value at log(shift_floor) instead of
% letting it run off to -Inf.  This also covers the case where min_cutoff sits
% only marginally below observed_min (gap positive but negligible), which a
% plain "shift = 0 when ~min_is_observed" rule would mishandle.
shift_floor = 1e-3 * max(spread, eps);
gap         = max(observed_min - min_cutoff, 0);
if gap < shift_floor
    shift = shift_floor - gap;
else
    shift = 0;
end

if do_log_transform
    P_trans = log(P - min_cutoff + shift);
    Q_trans = log(Q - min_cutoff + shift);
else
    P_trans = P;
    Q_trans = Q;
end

% --- Bandwidth extraction (observed samples) ---
% Bandwidths used for the KDE itself.  These MUST be obtained with the same
% Support option that will be used at evaluation time, because ksdensity
% interprets an explicit 'Bandwidth' on its internal working scale.
% These serve the point estimate d and the grid sizing below.  Bootstrap and
% permutation replicates do NOT reuse them -- score_pair re-estimates from each
% replicate's own data.
use_bounded = ~do_log_transform && ~min_is_observed;

if use_bounded
    [~, ~, bw_P] = ksdensity(P_trans, 'Support', [min_cutoff, Inf]);
    [~, ~, bw_Q] = ksdensity(Q_trans, 'Support', [min_cutoff, Inf]);
else
    [~, ~, bw_P] = ksdensity(P_trans);
    [~, ~, bw_Q] = ksdensity(Q_trans);
end

% --- PATCH (1): grid-extension bandwidth on the grid's own scale ----------
% bw_ext is added to grid coordinates, so it must share their units.  In the
% bounded-support branch bw_P/bw_Q live on ksdensity's internal
% log(x - min_cutoff) scale and would be meaningless added to raw data, so a
% raw-scale bandwidth is estimated separately for sizing purposes only.  In
% every other branch the KDE bandwidth is already on the grid scale (log space
% for the log path, raw space for the unbounded path) and is reused directly.
if use_bounded
    [~, ~, bw_grid_P] = ksdensity(P_trans);
    [~, ~, bw_grid_Q] = ksdensity(Q_trans);
else
    bw_grid_P = bw_P;
    bw_grid_Q = bw_Q;
end

% Extend grid 3 bandwidths beyond data range to capture KDE tail mass
bw_ext = 3 * max(bw_grid_P, bw_grid_Q);

if do_log_transform
    lo = min([P_trans; Q_trans]) - bw_ext;
    hi = max([P_trans; Q_trans]) + bw_ext;
else
    hi = max([P; Q]) + bw_ext;
    if use_bounded
        % --- PATCH (2): keep the first node strictly inside the support ---
        % Bounded-support KDE is singular at x = min_cutoff: the internal
        % transform gives log(0) = -Inf and the back-transform gives 0/0, i.e.
        % NaN.  That NaN was previously absorbed silently by max(p_pdf, tiny)
        % (MATLAB's max ignores NaN), quietly discarding the node.  The
        % estimated density tends to 0 as x -> min_cutoff, so nudging the first
        % node inside the support loses no mass and removes the singularity.
        lo = min_cutoff + 1e-6 * max(hi - min_cutoff, eps);
    else
        lo = min_cutoff;
    end
end
grid_eval = linspace(lo, hi, Ngrid);

if lo == hi
    d = 0; d_CI = []; d_SE = []; d_boot = []; noise_P = 0; noise_Q = 0;
    d_corr = 0; d_corr_CI = []; d_corr_SE = [];
    p_perm = NaN; d_perm = [];   % zero-range grid: nothing to test
    return;
end

%% ------------------ Main Computation ------------------
[d, div_val, ~] = compute_jsd_core(P_trans, Q_trans, grid_eval, do_log_transform, min_cutoff, min_is_observed, bw_P, bw_Q);

%% ------------------ Optimized Bootstrap Loop ------------------
d_CI = []; d_boot = []; d_SE = []; noise_P = 0; noise_Q = 0; d_corr = 0; d_corr_CI = []; d_corr_SE = [];
p_perm = NaN; d_perm = [];

if n_boot > 0
    div_boot   = zeros(n_boot,1);
    divP_noise = zeros(n_boot,1);
    divQ_noise = zeros(n_boot,1);

    n_degen_boot = 0;

    % Serial loop, drawing from the global stream in a fixed order. Do not
    % reorder these six randi calls or interpose others: the draw sequence IS
    % the reproducibility contract, and any change silently invalidates
    % bootstrap outputs saved from earlier runs at the same rng_seed.
    %
    % Bandwidth is re-estimated inside score_pair from each resample, not
    % inherited from the observed sample. That costs runtime and changes
    % d_CI/d_SE/noise_*/d_corr relative to versions of this function that
    % cached bandwidth, but it is what makes the noise floor an honest estimate
    % of the bias in d. See the Bandwidth note in ALGORITHM NOTES. The randi
    % sequence above is unaffected, so the resamples themselves are identical.
    for b = 1:n_boot
        % Resample indices
        P_star = P_trans(randi(nP, nP, 1));
        Q_star = Q_trans(randi(nQ, nQ, 1));
        P_s1   = P_trans(randi(nP, nP, 1));
        P_s2   = P_trans(randi(nP, nP, 1));
        Q_s1   = Q_trans(randi(nQ, nQ, 1));
        Q_s2   = Q_trans(randi(nQ, nQ, 1));

        [div_boot(b),   g1] = score_pair(P_star, Q_star, grid_eval, do_log_transform, min_cutoff, min_is_observed, use_bounded);
        [divP_noise(b), g2] = score_pair(P_s1,   P_s2,   grid_eval, do_log_transform, min_cutoff, min_is_observed, use_bounded);
        [divQ_noise(b), g3] = score_pair(Q_s1,   Q_s2,   grid_eval, do_log_transform, min_cutoff, min_is_observed, use_bounded);
        n_degen_boot = n_degen_boot + double(g1 || g2 || g3);
    end

    % --- Noise correction applied outside the loop ---
    avg_noise_div = 0.5 * (mean(divP_noise) + mean(divQ_noise));

    % Convert divergence arrays back to distances for outputs
    d_boot      = sqrt(max(0, div_boot));
    d_boot_corr = sqrt(max(0, div_boot - avg_noise_div));

    % Alpha Quantiles
    lo_q = ci_alpha/2;
    hi_q = 1 - ci_alpha/2;

    d_CI = quantile(d_boot, [lo_q, hi_q]);
    if n_boot > 1, d_SE = std(d_boot); else, d_SE = 0; end

    % Noise limits evaluated in Distance space at output
    noise_P = mean(sqrt(divP_noise));
    noise_Q = mean(sqrt(divQ_noise));

    % Final Divergence Noise Adjustment
    d_corr     = sqrt(max(0, div_val - avg_noise_div));
    d_corr_CI  = quantile(d_boot_corr, [lo_q, hi_q]);
    d_corr_SE  = std(d_boot_corr);

    if n_degen_boot > 0
        warning('jsd_kde:bootDegenerateReplicates', ...
            ['%d of %d bootstrap replicates drew a constant resample and ' ...
            'were scored by the closed-form limit; d_CI and the noise floor ' ...
            'are inflated accordingly.'], n_degen_boot, n_boot);
    end
end

%% ------------------ Label-Permutation Null ------------------
% Placed after the bootstrap loop on purpose: the bootstrap's draw sequence is
% a reproducibility contract (see the comment in that loop), so the permutation
% draws must not precede or interleave with it.  One randperm per replicate is
% this loop's own contract.
%
% Cost is roughly four ksdensity calls per replicate (two for bandwidth, two
% inside compute_jsd_core), so n_perm = 2000 is seconds-to-a-minute territory
% for the sample sizes this function is typically given.  The bootstrap costs
% three times that per replicate, so when both are on, n_boot usually wants to
% be the smaller of the two.
if n_perm > 0
    pool     = [P_trans; Q_trans];
    nTot     = nP + nQ;
    div_perm = zeros(n_perm, 1);
    n_degen  = 0;

    for b = 1:n_perm
        ix = randperm(nTot);
        A  = pool(ix(1:nP));
        B  = pool(ix(nP+1:nTot));

        % Transform, shift, cutoff and grid all held at their observed-data
        % values, so every replicate is scored by the same estimator as d.
        % Bandwidth is re-estimated per pseudo-group inside score_pair, at the
        % same group sizes as the observed statistic.
        [div_perm(b), g] = score_pair(A, B, grid_eval, do_log_transform, ...
            min_cutoff, min_is_observed, use_bounded);
        n_degen = n_degen + double(g);
    end

    d_perm = sqrt(max(0, div_perm));

    % Add-one estimator, compared in divergence space; the sqrt is monotone, so
    % the ordering and hence the p-value are identical either way.
    p_perm = (1 + sum(div_perm >= div_val)) / (1 + n_perm);

    if n_degen > 0
        warning('jsd_kde:permDegenerateReplicates', ...
            ['%d of %d permutation replicates produced a constant pseudo-' ...
            'group and were scored by the closed-form limit, which biases ' ...
            'p_perm upward (conservative).'], n_degen, n_perm);
    end
end
end

%% ============================================================
%  REPLICATE SCORER (shared by the bootstrap and permutation loops)
%% ============================================================

function [div, degenerate] = score_pair(A, B, grid_eval, do_log, min_cutoff, min_is_observed, use_bounded)
%SCORE_PAIR  Divergence for one resampled or permuted pair of samples.
%
%   Bandwidth is estimated from A and B themselves, so each replicate pays the
%   same bandwidth-selection cost as the observed statistic did at the same
%   sample sizes.  This is what lets the bootstrap noise floor stand in for the
%   bias in d; inheriting the observed sample's bandwidth understates it badly
%   when a sample is small.  See the Bandwidth note in ALGORITHM NOTES.
%
%   A resample or pseudo-group can be constant even when the observed sample is
%   not -- draw ten values from a tied-heavy sample and they can all coincide.
%   KDE carries no scale information there, so rather than let ksdensity
%   improvise a zero bandwidth this returns the same closed-form limits the
%   degenerate-input path uses: 0 for coincident point masses, 1 for mutually
%   singular ones.  DEGENERATE flags that case so callers can count and report
%   it, since those replicates sit at the top of the scale and make any
%   interval or p-value built from them conservative.

A_const = (max(A) - min(A)) <= 10*eps(max(abs(A)));
B_const = (max(B) - min(B)) <= 10*eps(max(abs(B)));

if A_const || B_const
    degenerate = true;
    if A_const && B_const && abs(A(1) - B(1)) <= 10*eps(max(abs([A(1), B(1)])))
        div = 0;   % identical point masses
    else
        div = 1;   % mutually singular; d = sqrt(1) = 1 on the bits basis
    end
    return;
end

degenerate = false;

% Bandwidth MUST be obtained with the same Support option used at evaluation
% time, because ksdensity interprets an explicit 'Bandwidth' on its internal
% working scale.  This mirrors the observed-sample extraction above.
if use_bounded
    [~, ~, bwA] = ksdensity(A, 'Support', [min_cutoff, Inf]);
    [~, ~, bwB] = ksdensity(B, 'Support', [min_cutoff, Inf]);
else
    [~, ~, bwA] = ksdensity(A);
    [~, ~, bwB] = ksdensity(B);
end

[~, div] = compute_jsd_core(A, B, grid_eval, do_log, min_cutoff, min_is_observed, bwA, bwB);
end

%% ============================================================
%  LIGHTWEIGHT CORE ESTIMATOR (Evaluated on Fixed Grid)
%% ============================================================

function [d, jsd_val, xgrid] = compute_jsd_core(P_t, Q_t, grid_eval, do_log, min_cutoff, min_is_observed, bw_P, bw_Q)
% NOTE: when do_log is true, xgrid (and hence the returned third output) is in
% LOG space, not on the original data scale.  See SCALE INVARIANCE above: the
% divergence is identical either way, so no back-transform is performed.
xgrid = grid_eval;

% JSD evaluated directly in the transformed space using invariance of
% f-divergences under monotone reparametrisation.  Avoiding the Jacobian
% back-transform keeps the grid uniform and removes the division by exp(y),
% both of which improve the accuracy of the trapezoidal quadrature.
if do_log
    p_pdf = ksdensity(P_t, xgrid, 'Function','pdf', 'Bandwidth', bw_P);
    q_pdf = ksdensity(Q_t, xgrid, 'Function','pdf', 'Bandwidth', bw_Q);
else
    if ~min_is_observed
        p_pdf = ksdensity(P_t, xgrid, 'Function','pdf', 'Support', [min_cutoff, Inf], 'Bandwidth', bw_P);
        q_pdf = ksdensity(Q_t, xgrid, 'Function','pdf', 'Support', [min_cutoff, Inf], 'Bandwidth', bw_Q);
    else
        p_pdf = ksdensity(P_t, xgrid, 'Function','pdf', 'Bandwidth', bw_P);
        q_pdf = ksdensity(Q_t, xgrid, 'Function','pdf', 'Bandwidth', bw_Q);
    end
end

% Explicit non-finite sanitisation.  max(x, tiny) already maps NaN to tiny
% because MATLAB's max ignores NaN, but it maps +Inf to +Inf, which would
% poison the quadrature.  Handle both cases the same way, before the floor.
p_pdf(~isfinite(p_pdf)) = 0;
q_pdf(~isfinite(q_pdf)) = 0;

% Subnormal-safe numerical floor
tiny = 1e-300;
p_pdf = max(p_pdf, tiny);
q_pdf = max(q_pdf, tiny);

% Normalize continuous densities safely
p_pdf = p_pdf / trapz(xgrid, p_pdf);
q_pdf = q_pdf / trapz(xgrid, q_pdf);

% Mixture distribution
m = 0.5 * (p_pdf + q_pdf);

% Logical guarding against 0/0 and 0*log(0)
klP_v = zeros(size(p_pdf));
klQ_v = zeros(size(q_pdf));

idxP = p_pdf > tiny;
idxQ = q_pdf > tiny;

klP_v(idxP) = p_pdf(idxP) .* log(p_pdf(idxP) ./ m(idxP));
klQ_v(idxQ) = q_pdf(idxQ) .* log(q_pdf(idxQ) ./ m(idxQ));

% Numerical integration on uniform grid
klP = trapz(xgrid, klP_v);
klQ = trapz(xgrid, klQ_v);

% Compute JSD in bits
jsd_val = max(0, 0.5 * (klP + klQ) / log(2));
d = sqrt(jsd_val);
end
