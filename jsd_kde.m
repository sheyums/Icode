function [d, d_CI, d_SE, d_boot, noise_P, noise_Q, d_corr, d_corr_CI, d_corr_SE] = jsd_kde(P, Q, Ngrid, min_cutoff, do_log_transform, n_boot, ci_alpha)
% JSD_KDE  Jensen-Shannon distance between two univariate distributions
%          estimated via kernel density estimation (KDE).
%
% SYNTAX
%   d = jsd_kde(P, Q)
%   [d, d_CI, d_SE, d_boot, noise_P, noise_Q, d_corr, d_corr_CI, d_corr_SE] = ...
%       jsd_kde(P, Q, Ngrid, min_cutoff, do_log_transform, n_boot, ci_alpha)
%
% SCALE INVARIANCE
%   JS divergence is exactly invariant under a smooth invertible
%   reparametrisation y = g(x): the Jacobian cancels inside the log because
%   the mixture transforms identically to its components,
%       p_Y/m_Y = (p_X/g')/(m_X/g') = p_X/m_X,
%   and the remaining g' cancels against dy = g' dx.  The estimator therefore
%   evaluates the integrals directly in the transformed space on a uniform
%   grid rather than back-transforming to the original scale.
%
% PATCH NOTES (relative to the previous revision)
%   (1) BANDWIDTH UNITS.  Bounded-support ksdensity works internally on
%       log(x - min_cutoff) and returns a bandwidth on THAT scale, so it
%       cannot be added to raw data values when sizing the evaluation grid.
%       A separate raw-scale bandwidth is now obtained for the grid extension
%       while the KDE itself keeps the bounded-support bandwidth.
%   (2) BOUNDARY OFFSET.  With bounded support the KDE is singular at exactly
%       x = min_cutoff (internal log(0) -> -Inf, then 0/0).  The first grid
%       node is now placed just inside the support.
%   (3) CONDITIONAL SHIFT.  The log-transform offset is applied only to the
%       extent needed to keep the smallest observation a sensible distance
%       above min_cutoff.  When min_cutoff already sits comfortably below the
%       data the shift is exactly zero, leaving log(x - min_cutoff) unperturbed.

%% ------------------ Defaults & Setup ------------------
arguments
    P           {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeFinite}
    Q           {mustBeNumeric, mustBeVector, mustBeNonempty, mustBeFinite}
    Ngrid       (1,1) {mustBeInteger, mustBePositive}                       = 512
    min_cutoff  (1,1) double                                                = NaN
    do_log_transform                                                        = []
    n_boot      (1,1) {mustBeInteger, mustBeNonnegative}                    = 0
    ci_alpha    (1,1) {mustBeInRange(ci_alpha, 0, 1, "exclusive")}          = 0.05
end

P = P(:); Q = Q(:);
nP = numel(P); nQ = numel(Q);

% Resolve data-dependent default
if isnan(min_cutoff), min_cutoff = min([P; Q]); end

if any(P < min_cutoff) || any(Q < min_cutoff)
    error('jsd_kde:belowCutoff', 'Values fall below min_cutoff.');
end

observed_min = min([P;Q]);
min_is_observed = abs(min_cutoff - observed_min) < 10*eps(observed_min);

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
