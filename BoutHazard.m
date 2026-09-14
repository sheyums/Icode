function H = BoutHazard(eventseries, xmin, options)
%BOUTHAZARD Life-table hazard with confidence bands, for left-truncated bouts.
%
%   H = BOUTHAZARD(eventseries, xmin, SamplingInterval=dt, ...)
%   H = BOUTHAZARD(..., SurvivalHandle=R.Fits(i).SurvivalHandle)
%   H = BOUTHAZARD(..., Bootstrap=999, Plot=true)
%
%   WHAT THIS IS FOR, AND WHAT IT IS NOT FOR. It does not rank models --
%   AICc, BIC, the G test and the bootstrap LRT do that, and the hazard is
%   a transform of the same fit, h = f/S, so it adds nothing to a
%   comparison. What it does is decide whether a whole FAMILY can work at
%   all, which is upstream of choosing within one:
%
%     - A hyperexponential's hazard is strictly decreasing at EVERY order:
%       a mixture of exponentials is a sum of decreasing terms. So is every
%       other monotone family. An exponentiated Weibull allows at most ONE
%       turning point.
%     - Data whose hazard falls, rises, then falls again are therefore
%       outside all of them BY CONSTRUCTION rather than by evidence, and no
%       information criterion will tell you that -- it can only rank the
%       candidates you thought of.
%     - Pooling individuals cannot rescue it either. h_mix(t) = sum_i
%       w_i(t) h_i(t) with the weights shifting toward the longer-lived
%       components, so mixing makes a hazard fall FASTER; it cannot create
%       a rise. A hump is therefore not evidence of heterogeneity.
%
%   WHY THE BANDS ARE THE POINT. A hump read off an unsmoothed plot is not
%   a finding. The risk set thins as t grows, so the tail of any hazard
%   estimate is noisy exactly where humps tend to appear, and a rise of
%   50% over a decade can be nothing. The question this answers is whether
%   a non-increasing hazard fits INSIDE the bands. If it does, the
%   family-exclusion argument above is not available and the monotone
%   families were never ruled out.
%
%   TRUNCATION IS FREE HERE, which is why the hazard is the right statistic
%   for this data. h(t) = P(end in (t, t+dt] | survived past t) is already
%   CONDITIONAL, so restricting to bouts >= xmin leaves it unbiased above
%   xmin -- no S(xmin) normalizer, none of the extrapolation that makes the
%   untruncated weight w treacherous. The estimator is the actuarial one,
%
%       h_j = d_j / n_j,   n_j = #{t > e_j},  d_j = #{e_j < t <= e_{j+1}}
%
%   converted to a RATE per unit time as
%
%       rate_j = -log(1 - h_j) / (e_{j+1} - e_j)
%
%   and NOT as h_j / width, which is wrong for wide bins and silently so.
%   Under a constant hazard lam across a bin of width W the interval
%   probability is p = 1 - exp(-lam W), so p/W tends to lam only as W -> 0;
%   with log-spaced bins the tail bins are decades wide, p -> 1, and p/W
%   collapses toward 1/W regardless of the true rate. On exponential data
%   that error read as a hazard falling by a factor of three across bins
%   where the truth is flat -- a spurious DECREASE, which is exactly the
%   shape this function exists to distinguish from a real one. The log
%   form inverts the relation exactly and is the standard actuarial
%   conversion (constant hazard within the interval).
%
%   CONFIDENCE BANDS are Wilson score intervals on the binomial d_j/n_j
%   (Wilson 1927), not Wald. Wald intervals on a proportion collapse to
%   zero width as d_j -> 0 and can cover below zero, which is precisely the
%   regime of a hazard's tail; Brown, Cai & DasGupta (2001) is the standard
%   reference for why the textbook interval should not be used here. The
%   bands are pointwise, NOT simultaneous: reading "some bin rises" off 24
%   of them invites a multiplicity error, which is what RiseRatio and its
%   bootstrap interval are for.
%
%   BINS default to log-spaced, since bout durations span decades. A bin
%   left with fewer than MinAtRisk survivors is dropped rather than plotted
%   as noise, and so is a SATURATED bin in which every survivor ended --
%   its rate is unbounded, not large, and the final bin is always saturated
%   because it runs to max(t). H.nDropped and H.nSaturated record both.
%
%   OPTIONS
%     SamplingInterval  required, same units as xmin, as everywhere here
%     BinEdges          explicit edges; overrides NumBins/Spacing
%     NumBins           default 24
%     Spacing           "log" (default) or "linear"
%     MinAtRisk         default 10; bins thinner than this are dropped
%     Alpha             default 0.05
%     SurvivalHandle    a fitted TRUNCATED survival S(t)/S(xmin), e.g. from
%                       CompareBoutModels. Its hazard is computed on the
%                       SAME bins by the same formula, so the comparison is
%                       like for like rather than a continuous curve
%                       against a binned estimate.
%     Bootstrap         B resamples for the RiseRatio interval; 0 skips
%     NullSurvivalHandle  a fitted MONOTONE model's truncated survival, to
%                       test the rise against. Without it there is no null
%                       and NonMonotone falls back to the band-overlap rule
%     NullReplicates    draws from that null; nothing is refitted, so this
%                       is cheap -- 999 costs seconds
%     RandomSeed        pinned to 'twister', as elsewhere here
%     Plot              hazard with bands, log-log
%
%   OUTPUT
%     H.BinEdges, H.BinCenters, H.BinWidth
%     H.AtRisk, H.Events          n_j and d_j
%     H.Hazard, H.Lower, H.Upper  per unit time
%     H.ModelHazard               same bins, if SurvivalHandle was given
%     H.TroughIndex, H.PeakIndex  minimum (excluding the final bin, which
%                                 has nothing after it), and the maximum
%                                 after it
%     H.RiseRatio                 h_peak / h_trough. >= 1 for ANY curve, so
%                                 large is only meaningful against a null
%     H.RiseCI                    bootstrap percentile interval -- DESCRIPTIVE,
%                                 not a test; see the note in the code
%     H.RiseNullP                 p-value against NullSurvivalHandle, the
%                                 only calibrated statement here
%     H.RiseNullQuantiles         median and 95th percentile of the null
%     H.RiseDisjoint              peak's lower band above trough's upper
%     H.NonMonotone               RiseNullP < Alpha when a null was given --
%                                 the only calibrated form. Falls back to
%                                 RiseDisjoint when it was not. It is NOT
%                                 "RiseCI excludes 1": against 50 datasets
%                                 with a strictly decreasing hazard that
%                                 fired 20% of the time at a nominal 5%.
%     H.n, H.xmin, H.SamplingInterval, H.nDropped, H.nSaturated
%
%   See also COMPAREBOUTMODELS, FITHYPERERLANGMLE, FITHYPEREXPONENTIALMLE.

arguments
    eventseries double {mustBeReal}
    xmin (1,1) double {mustBePositive}
    options.SamplingInterval (1,1) double = NaN
    options.BinEdges double = []
    options.NumBins (1,1) double {mustBeInteger,mustBePositive} = 24
    options.Spacing (1,1) string {mustBeMember(options.Spacing,["log","linear"])} = "log"
    options.MinAtRisk (1,1) double {mustBeInteger,mustBePositive} = 10
    options.Alpha (1,1) double {mustBePositive} = 0.05
    options.SurvivalHandle = []
    options.NullSurvivalHandle = []
    options.NullReplicates (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.RandomSeed = []
    options.Plot (1,1) logical = false
end

if isnan(options.SamplingInterval)
    error('BoutHazard:SamplingIntervalRequired', ...
        ['SamplingInterval is required, in the SAME UNITS as xmin. It sets ' ...
         'the grid the bin edges are snapped to.']);
elseif ~(options.SamplingInterval > 0)
    error('BoutHazard:InvalidSamplingInterval', ...
        'SamplingInterval must be a positive scalar.');
end
dt = options.SamplingInterval;

t = eventseries(:);
t = t(isfinite(t) & t >= xmin);
n = numel(t);
if n < 2
    error('BoutHazard:TooFewData', ...
        'Only %d bouts at or above xmin=%g; nothing to estimate.', n, xmin);
end

% ------------------------------------------------------------------- bins
edges = buildEdges(t, xmin, dt, options);
if numel(edges) < 2
    error('BoutHazard:DegenerateBins', ...
        ['Bin edges collapsed to %d distinct grid points. The data span ' ...
         'less than one SamplingInterval, or NumBins is far too large.'], ...
        numel(edges));
end

nb = numel(edges) - 1;
atRisk = zeros(1, nb);
events = zeros(1, nb);
for j = 1:nb
    atRisk(j) = sum(t > edges(j));
    events(j) = sum(t > edges(j) & t <= edges(j+1));
end

width = diff(edges);

% A SATURATED bin -- every survivor ended inside it, d_j = n_j -- has an
% unbounded rate, not a large one: -log(1-p)/W is infinite at p = 1. The
% final bin is saturated BY CONSTRUCTION, since it runs to max(t) and
% nobody can survive past the longest observed bout, which is why the
% actuarial convention leaves the last interval open-ended and does not
% quote a hazard for it. Excluded rather than reported as Inf: left in, it
% became the peak of every bootstrap replicate, so every RiseRatio came
% back Inf and the interval was NaN on data with an obvious rise.
saturated = events >= atRisk;
keep = atRisk >= options.MinAtRisk & ~saturated;
nDropped = sum(~keep);
nSaturated = sum(saturated);

edgesK = edges;
atRisk = atRisk(keep); events = events(keep); width = width(keep);
lefts = edgesK(1:end-1); rights = edgesK(2:end);
lefts = lefts(keep); rights = rights(keep);
if isempty(atRisk)
    error('BoutHazard:NoUsableBins', ...
        ['Every bin held fewer than MinAtRisk=%d survivors. Use fewer bins, ' ...
         'or lower MinAtRisk and read the tail with suspicion.'], ...
        options.MinAtRisk);
end

p = events ./ atRisk;
[plo, phi] = wilson(p, atRisk, options.Alpha);

% Interval probability -> rate. Monotone in p, so the Wilson endpoints map
% straight through and the band needs no separate derivation. p = 1 means
% every survivor ended inside the bin, where the rate is unbounded: it is
% reported as Inf rather than clamped, so it cannot masquerade as a finite
% estimate, and it is excluded from the rise statistics below.
toRate = @(pp, ww) -log1p(-min(pp, 1)) ./ ww;

H = struct();
H.BinEdges    = edgesK;
H.BinCenters  = sqrt(max(lefts, dt/2) .* rights);   % geometric, for log bins
H.BinWidth    = width;
H.AtRisk      = atRisk;
H.Events      = events;
H.Hazard      = toRate(p, width);
H.Lower       = toRate(plo, width);
H.Upper       = toRate(phi, width);
H.nDropped    = nDropped;
H.nSaturated  = nSaturated;      % of those, how many had d_j = n_j
H.n           = n;
H.xmin        = xmin;
H.SamplingInterval = dt;
H.Alpha       = options.Alpha;

% ------------------------------------------------- the fitted model, if any
H.ModelHazard = [];
if ~isempty(options.SurvivalHandle)
    S = options.SurvivalHandle;
    sl = S(lefts); sr = S(rights);
    pm = zeros(size(sl));
    ok = sl > 0;
    pm(ok) = (sl(ok) - sr(ok)) ./ sl(ok);      % same estimand as d_j/n_j
    pm(~ok) = NaN;
    H.ModelHazard = toRate(pm, width);         % ...and the same conversion
end

% ------------------------------------------------------------- the rise
% The trough is searched among bins that HAVE a successor. Without that,
% the global minimum can land on the last kept bin -- a thin tail bin with
% a wide band -- leaving trough = peak and RiseRatio = 1, a false NEGATIVE
% produced by one noisy bin. Seen on real bouts: at 16 bins the statistic
% found a rise of 1.39, at 20 bins the minimum moved to the final bin
% (35 at risk) and the same data reported no rise at all.
hh = H.Hazard; hh(~isfinite(hh)) = Inf;
if numel(hh) >= 2
    [~, iT] = min(hh(1:end-1));
else
    iT = 1;
end
H.TroughIndex = iT;
if iT < numel(hh)
    tail = hh(iT+1:end); tail(~isfinite(tail)) = -Inf;
    [~, rel] = max(tail);
    iP = iT + rel;
else
    iP = iT;
end
if ~isfinite(H.Hazard(iP)), iP = iT; end
H.PeakIndex = iP;
if H.Hazard(iT) > 0
    H.RiseRatio = H.Hazard(iP) / H.Hazard(iT);
else
    H.RiseRatio = Inf;
end
H.RiseDisjoint = iP > iT && H.Lower(iP) > H.Upper(iT);

% Bootstrap the ratio over BOUTS, recomputing on the same bins.
%
% READ THIS AS A DESCRIPTIVE INTERVAL, NOT A TEST. RiseRatio is a maximum
% taken after a minimum, so it is >= 1 for any curve whatsoever, and the
% bootstrap resamples the DATA rather than drawing from a no-rise null.
% Its lower limit therefore exceeds 1 for almost any noisy hazard, monotone
% or not: on 50 datasets simulated from a strictly DECREASING hazard,
% RiseCI(1) > 1 held in 9 of 50. It says how precisely the observed rise is
% measured; it does not say the rise is real.
%
% For that, supply NullSurvivalHandle -- a fitted MONOTONE model -- and
% RiseNullP below compares the observed ratio against ratios from data
% simulated under it. That is the comparison with a null in it.
H.RiseCI = [NaN NaN];
if options.Bootstrap > 0
    if ~isempty(options.RandomSeed)
        rng(options.RandomSeed, 'twister');   % 'twister' is not decoration;
    end                                        % see CLAUDE.md
    B = options.Bootstrap;
    rb = nan(1, B);
    for b = 1:B
        tb = t(randi(n, n, 1));
        hb = binHazard(tb, edgesK, keep, width);
        if isempty(hb) || ~any(isfinite(hb)), continue; end
        rb(b) = riseOf(hb);
    end
    rb = rb(isfinite(rb));
    if ~isempty(rb)
        H.RiseCI = quantilePct(rb, [100*options.Alpha/2, 100*(1-options.Alpha/2)]);
    end
end

% NonMonotone is RiseDisjoint ALONE. It was once RiseDisjoint OR
% RiseCI(1) > 1, and that second clause is not a test: see the RiseCI
% comment above. Measured against 50 datasets simulated from a strictly
% decreasing hazard, the pair fired 20% of the time while RiseDisjoint
% alone fired 6%, against a nominal 5%.
% ------------------------------------------------- the test WITH a null
% RiseRatio on its own has no null: it is >= 1 for every curve, so there
% is nothing for it to be large RELATIVE TO. This supplies one. Simulate
% datasets of the same size from a fitted MONOTONE model -- a
% hyperexponential, whose hazard is strictly decreasing at any order --
% recompute the same statistic on the same bins, and ask how often the
% null reaches the observed value:
%
%     p = (1 + #{ratio_null >= ratio_observed}) / (B + 1)
%
% Nothing is refitted, so this is cheap: the replicates are draws from a
% survival curve, not maximum-likelihood fits. It is the same logic as
% HYPEREXPONENTIALLRT's bootstrap, for a different statistic.
H.RiseNullP = NaN;
H.RiseNullQuantiles = [NaN NaN];
if ~isempty(options.NullSurvivalHandle) && options.NullReplicates > 0
    if ~isempty(options.RandomSeed)
        rng(options.RandomSeed + 1, 'twister');   % not the resampling stream
    end
    Bn = options.NullReplicates;
    rn = nan(1, Bn);
    % The grid holds the ATTAINABLE durations, xmin upward. The engine's
    % truncated survival is 1 at xmin-dt, NOT at xmin, so the mass landing
    % on g_k is S(g_k - dt) - S(g_k). Differencing S(grid) instead put
    % every atom one step early: measured against the engine's own pmf,
    % E[T] came out 0.99 s low and P(T=xmin) equalled the engine's
    % P(xmin) + P(xmin+dt).
    grid_ = (xmin) : dt : max(edgesK(end), max(t));
    Sprev = options.NullSurvivalHandle(grid_ - dt);
    Snow  = options.NullSurvivalHandle(grid_);
    Sprev = min(max(Sprev(:).', 0), 1);
    Snow  = min(max(Snow(:).',  0), 1);
    Sprev(1) = 1;                    % S(xmin - dt) = 1 by construction
    pmf = max(Sprev - Snow, 0);
    pmf(end) = pmf(end) + max(Snow(end), 0);   % residual tail as one atom
    tot = sum(pmf);
    if tot > 0
        cdf_ = cumsum(pmf / tot);
        for b = 1:Bn
            u = rand(n, 1);
            idx = arrayfun(@(uu) find(cdf_ >= uu, 1), u);
            tb = grid_(idx).';
            % EACH REPLICATE GETS ITS OWN BINS. Reusing the observed data's
            % edges and keep mask makes the null and the observed statistic
            % non-exchangeable: the observed ratio was computed after its
            % own bin and at-risk selection, so holding those fixed for the
            % null removes a source of variation the observed value had.
            % Measured on 999 draws from a fitted K=2: fixing the bins
            % raised p from 0.019 to 0.047 at 16 bins and from 0.057 to
            % 0.185 at 20 -- a factor of 2.5 to 3, and the whole of the
            % conservatism. The sampler shift above changed nothing.
            hb = hazardOf(tb, xmin, dt, options);
            rn(b) = riseOf(hb);
        end
        rn = rn(isfinite(rn));
        if ~isempty(rn)
            H.RiseNullP = (1 + sum(rn >= H.RiseRatio)) / (numel(rn) + 1);
            H.RiseNullQuantiles = quantilePct(rn, [50, 95]);
        end
    end
end

% NonMonotone is RiseDisjoint ALONE when no null was supplied. It was once
% RiseDisjoint OR RiseCI(1) > 1, and that second clause is not a test: see
% the RiseCI comment above. Measured against 50 datasets simulated from a
% strictly decreasing hazard, the pair fired 20% of the time while
% RiseDisjoint alone fired 6%, against a nominal 5%. With a null supplied,
% RiseNullP decides it, which is the only form of this that has a
% calibrated error rate.
if isfinite(H.RiseNullP)
    H.NonMonotone = H.RiseNullP < options.Alpha;
else
    H.NonMonotone = H.RiseDisjoint;
end

if options.Plot
    H.Figure = plotHazard(H);
else
    H.Figure = [];
end
end

% ------------------------------------------------------------------------
function edges = buildEdges(t, xmin, dt, options)
%BUILDEDGES The bin edges for one sample. Shared by the observed data and
%by every null replicate, so a replicate is binned the way the data were.
if ~isempty(options.BinEdges)
    edges = sort(options.BinEdges(:).');
else
    lo = xmin - dt;                    % so bouts exactly at xmin are events
    hi = max(t);
    if options.Spacing == "log"
        % log spacing over the RETAINED range. lo can be 0 when xmin == dt,
        % so the log grid starts at the first positive edge and lo is
        % prepended, rather than taking log(0).
        first = max(lo, dt/2);
        e = exp(linspace(log(first), log(hi), options.NumBins + 1));
        edges = [lo, e(2:end)];
    else
        edges = linspace(lo, hi, options.NumBins + 1);
    end
end
% Snap to the sampling grid and drop duplicates: two edges inside one grid
% step describe a bin no observation can fall in, which would read as a
% hazard of exactly zero rather than as an empty bin.
edges = unique(round(edges / dt) * dt);
end

% ------------------------------------------------------------------------
function h = hazardOf(t, xmin, dt, options)
%HAZARDOF The whole estimator for one sample -- its own edges, its own
%at-risk and saturation screening, its own rate conversion. This is what a
%null replicate must go through for its statistic to be comparable with the
%observed one.
h = [];
edges = buildEdges(t, xmin, dt, options);
if numel(edges) < 2, return; end
nb = numel(edges) - 1;
aR = zeros(1, nb); ev = zeros(1, nb);
for j = 1:nb
    aR(j) = sum(t > edges(j));
    ev(j) = sum(t > edges(j) & t <= edges(j+1));
end
kp = aR >= options.MinAtRisk & ~(ev >= aR);
if ~any(kp), return; end
w = diff(edges);
h = -log1p(-min(ev(kp) ./ aR(kp), 1)) ./ w(kp);
end

% ------------------------------------------------------------------------
function r = riseOf(h)
%RISEOF max-after-min, the same statistic for the data, the bootstrap and
%the null -- so they are comparable by construction rather than by care.
%The trough excludes the final bin, which has nothing after it.
r = NaN;
h = h(isfinite(h));
if numel(h) < 2, return; end
[~, iT] = min(h(1:end-1));
hp = max(h(iT+1:end));
if h(iT) > 0, r = hp / h(iT); end
end

% ------------------------------------------------------------------------
function [lo, hi] = wilson(p, n, alpha)
%WILSON Score interval for a binomial proportion.
% Wald (p +- z sqrt(p(1-p)/n)) is wrong in exactly the regime a hazard tail
% lives in: at d = 0 it has zero width, and near the boundary it covers
% below 0. The score interval stays inside [0,1] and keeps its coverage at
% small counts (Brown, Cai & DasGupta 2001).
z = sqrt(2) * erfinv(1 - alpha);        % no norminv: no Statistics Toolbox
den = 1 + z.^2 ./ n;
ctr = (p + z.^2 ./ (2*n)) ./ den;
half = (z ./ den) .* sqrt(p .* (1-p) ./ n + z.^2 ./ (4 * n.^2));
lo = max(ctr - half, 0);
hi = min(ctr + half, 1);
end

% ------------------------------------------------------------------------
function h = binHazard(t, edges, keep, width)
%BINHAZARD The estimator again, for one bootstrap replicate, on fixed bins.
nb = numel(edges) - 1;
p = zeros(1, nb);
for j = 1:nb
    nr = sum(t > edges(j));
    if nr > 0
        p(j) = sum(t > edges(j) & t <= edges(j+1)) / nr;
    else
        p(j) = NaN;
    end
end
h = -log1p(-min(p(keep), 1)) ./ width;   % same conversion as the point estimate
end

% ------------------------------------------------------------------------
function q = quantilePct(x, pct)
%QUANTILEPCT Percentiles without the Statistics Toolbox.
x = sort(x(:));
m = numel(x);
q = zeros(1, numel(pct));
for i = 1:numel(pct)
    r = pct(i)/100 * (m - 1) + 1;
    lo = floor(r); hi = ceil(r);
    if lo == hi
        q(i) = x(lo);
    else
        q(i) = x(lo) + (r - lo) * (x(hi) - x(lo));
    end
end
end

% ------------------------------------------------------------------------
function f = plotHazard(H)
f = figure('Color', 'w');
% A new figure's THEME is not predictable in MATLAB R2026a: two figures
% created moments apart in one session came out light and dark, and an
% earlier session got the reverse. 'Color','w' does not settle it -- it
% whitens the figure margin while the axes stay dark (axes Color 0.07,
% XColor 0.85), so a dark data curve on a dark axes is invisible. Forced
% before the axes exist, so children inherit it. theme() is absent in
% Octave and pre-R2025a MATLAB, hence the try.
try, theme(f, 'light'); catch, end %#ok<TRYNC>
ax = axes(f); hold(ax, 'on');
xv = H.BinCenters;
good = isfinite(H.Hazard) & H.Hazard > 0;
patch(ax, [xv(good), fliplr(xv(good))], ...
          [max(H.Lower(good), realmin), fliplr(max(H.Upper(good), realmin))], ...
      [0.85 0.88 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
plot(ax, xv(good), H.Hazard(good), 'o-', 'Color', [0.15 0.25 0.55], ...
     'MarkerFaceColor', [0.15 0.25 0.55], 'MarkerSize', 4, 'LineWidth', 1.2);
if ~isempty(H.ModelHazard)
    gm = isfinite(H.ModelHazard) & H.ModelHazard > 0;
    plot(ax, xv(gm), H.ModelHazard(gm), '-', 'Color', [0.80 0.33 0.15], 'LineWidth', 1.6);
end
set(ax, 'XScale', 'log', 'YScale', 'log', 'Box', 'off');
xlabel(ax, 'bout duration'); ylabel(ax, 'hazard (per unit time)');
ttl = sprintf('n = %d, xmin = %g', H.n, H.xmin);
if H.NonMonotone
    ttl = [ttl, sprintf('  --  RISE x%.2f, bands disjoint', H.RiseRatio)];
else
    ttl = [ttl, sprintf('  --  rise x%.2f, not resolved', H.RiseRatio)];
end
title(ax, ttl);
end
