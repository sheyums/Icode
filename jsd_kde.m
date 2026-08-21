function [d, d_CI, d_SE, d_boot, noise_P, noise_Q, d_corr, d_corr_CI, d_corr_SE] = jsd_kde(P, Q, Ngrid, min_cutoff, do_log_transform, n_boot, ci_alpha, rng_seed)
% JSD_KDE  Jensen-Shannon distance between two univariate distributions
%          estimated via kernel density estimation (KDE).
%
% SYNTAX
%   d = jsd_kde(P, Q)
%   [d, d_CI, d_SE, d_boot, noise_P, noise_Q, d_corr, d_corr_CI, d_corr_SE] = ...
%       jsd_kde(P, Q, Ngrid, min_cutoff, do_log_transform, n_boot, ci_alpha, rng_seed)
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
% INPUTS
%   P, Q              Numeric column (or row) vectors of observations.
%                     All values must be >= min_cutoff.  Neither may be empty.
%   Ngrid             Number of evaluation grid points (default: 512).
%   min_cutoff        Hard lower bound of the support.  Defaults to
%                     min([P; Q]).  Set to a known physical minimum (e.g. 0 for
%                     positive-definite quantities) to enable bounded-support
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
%                     For publication-quality CI and noise correction,
%                     1000-2000 is a reasonable choice.
%   ci_alpha          Nominal coverage level, alpha in (0,1)
%                     (default: 0.05 => 95% CI).
%   rng_seed          Nonnegative integer seed for the bootstrap, or [] to draw
%                     from the current stream (default: []).  See
%                     REPRODUCIBILITY below.
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
%   d_corr_CI         Bootstrap CI on d_corr; [] if n_boot = 0.
%   d_corr_SE         Bootstrap standard error of d_corr; [] if n_boot = 0.
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
%   Bandwidth caching.  bw_P and bw_Q are extracted via ksdensity once before
%   any bootstrap iteration and passed explicitly to all internal ksdensity
%   calls, eliminating O(6 * n_boot) redundant bandwidth estimations that would
%   otherwise dominate runtime.  The bootstrap loop is a plain for-loop;
%   replacing it with parfor (Parallel Computing Toolbox) requires no further
%   changes.
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
%   The global stream is seeded on entry and its previous state is restored on
%   exit via onCleanup, including on error or Ctrl-C, so calling jsd_kde never
%   perturbs the caller's random stream.  Omit rng_seed (default []) to draw
%   from the current stream and leave it advanced.
%
%   Only d is deterministic without a seed.  d_CI, d_SE, d_boot, d_corr,
%   d_corr_CI and d_corr_SE all vary run to run, so results intended to be
%   reproducible from saved output should always be generated with a seed.
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

%% ------------------ Defaults & Setup ------------------
arguments
    P           {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeFinite}
    Q           {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeFinite}
    Ngrid       (1,1) {mustBeInteger, mustBePositive}                       = 512
    min_cutoff  (1,1) double                                                = NaN
    do_log_transform                                                        = []
    n_boot      (1,1) {mustBeInteger, mustBeNonnegative}                    = 0
    ci_alpha    (1,1) {mustBeInRange(ci_alpha, 0, 1, "exclusive")}          = 0.05
    rng_seed                                                                = []
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
    cleanup_rng  = onCleanup(@() rng(rng_state_in));  %#ok<NASGU>
    rng(rng_seed, 'twister');
end

% Resolve data-dependent default
if isnan(min_cutoff), min_cutoff = min([P; Q]); end

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
    return;
end

%% ------------------ Auto-log decision ------------------
% Auto-log functions when min_cutoff is at its default (observed minimum).
if isempty(do_log_transform)
    skew_thresh = 0.8 + 10/sqrt(nP + nQ);
    do_log_transform = (all(P >= min_cutoff) && all(Q >= min_cutoff)) && ...
        ((skewness(P) > skew_thresh) || (skewness(Q) > skew_thresh));
    if do_log_transform
        warning('jsd_kde:logTransform', ...
            'Log-transform applied automatically (skewness P=%.2f, Q=%.2f exceeded threshold %.2f).', ...
            skewness(P), skewness(Q), skew_thresh);
    end
elseif do_log_transform
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

% --- Bandwidth extraction (cached for bootstrap) ---
% Bandwidths used for the KDE itself.  These MUST be obtained with the same
% Support option that will be used at evaluation time, because ksdensity
% interprets an explicit 'Bandwidth' on its internal working scale.
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
    d_corr = 0; d_corr_CI = []; d_corr_SE = []; return;
end

%% ------------------ Main Computation ------------------
[d, div_val, ~] = compute_jsd_core(P_trans, Q_trans, grid_eval, do_log_transform, min_cutoff, min_is_observed, bw_P, bw_Q);

%% ------------------ Optimized Bootstrap Loop ------------------
d_CI = []; d_boot = []; d_SE = []; noise_P = 0; noise_Q = 0; d_corr = 0; d_corr_CI = []; d_corr_SE = [];

if n_boot > 0
    div_boot   = zeros(n_boot,1);
    divP_noise = zeros(n_boot,1);
    divQ_noise = zeros(n_boot,1);

    for b = 1:n_boot
        % Resample indices
        P_star = P_trans(randi(nP, nP, 1));
        Q_star = Q_trans(randi(nQ, nQ, 1));
        P_s1   = P_trans(randi(nP, nP, 1));
        P_s2   = P_trans(randi(nP, nP, 1));
        Q_s1   = Q_trans(randi(nQ, nQ, 1));
        Q_s2   = Q_trans(randi(nQ, nQ, 1));

        [~, div_boot(b)]   = compute_jsd_core(P_star, Q_star, grid_eval, do_log_transform, min_cutoff, min_is_observed, bw_P, bw_Q);
        [~, divP_noise(b)] = compute_jsd_core(P_s1,   P_s2,   grid_eval, do_log_transform, min_cutoff, min_is_observed, bw_P, bw_P);
        [~, divQ_noise(b)] = compute_jsd_core(Q_s1,   Q_s2,   grid_eval, do_log_transform, min_cutoff, min_is_observed, bw_Q, bw_Q);
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
    d_SE = std(d_boot);

    % Noise limits evaluated in Distance space at output
    noise_P = mean(sqrt(divP_noise));
    noise_Q = mean(sqrt(divQ_noise));

    % Final Divergence Noise Adjustment
    d_corr     = sqrt(max(0, div_val - avg_noise_div));
    d_corr_CI  = quantile(d_boot_corr, [lo_q, hi_q]);
    d_corr_SE  = std(d_boot_corr);
end
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
