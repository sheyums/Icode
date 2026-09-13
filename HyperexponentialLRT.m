function L = HyperexponentialLRT(eventseries, xmin, K0, K1, options)
%HYPEREXPONENTIALLRT Parametric bootstrap likelihood-ratio test for the
%order of a left-truncated hyperexponential mixture.
%
%   L = HYPEREXPONENTIALLRT(data, xmin, K0, K1, SamplingInterval=dt, ...)
%
%   Tests H0: K0 components against H1: K1 components (K0 < K1, so the
%   models are nested) by simulating the null distribution of the
%   likelihood-ratio statistic rather than assuming one.
%
%   WHY NOT CHI-SQUARE. For nested models the usual result is
%   LR = 2*(logL1 - logL0) ~ chi2 with df = the difference in free
%   parameters, here 2*(K1-K0). That result requires the null to be an
%   INTERIOR point of the alternative's parameter space with a
%   non-singular information matrix, and a mixture violates both. A K0
%   mixture sits inside a K1 mixture in two different degenerate ways --
%   an extra component with zero weight (a boundary point: weights cannot
%   go negative), or two components sharing a time constant (a whole flat
%   ridge, since mass slides between the twins with no change in
%   likelihood). Neither is an interior maximum, the information matrix is
%   singular along the ridge, and the chi-square approximation does not hold.
%   The true null distribution instead carries an ATOM AT ZERO, from the
%   replicates in which the extra component collapses.
%
%   The DIRECTION of the chi-square error is not guaranteed a priori, and
%   depends on the family: for normal mixtures with unrestricted variances
%   the LRT can diverge without bound (Hartigan 1985), the opposite
%   direction entirely. That is exactly why the null has to be simulated
%   rather than argued.
%
%   For THIS family, measured: on truncated exponential mixtures at n=800
%   the simulated null put 17-24% of its mass exactly at zero and had a
%   95th percentile of 4.0-4.7, against chi2(2)'s 5.99. The null is
%   therefore stochastically SMALLER, so chi-square is CONSERVATIVE here:
%   its critical value is too high and it under-rejects a real extra
%   component. Do not carry that direction over to another family, another
%   n, or another truncation -- read L.LRNull for the null actually
%   obtained in your case. See Hartigan (1985), McLachlan (1987) Appl
%   Statist 36:318-324, Lindsay (1995), and McLachlan & Peel (2000) Finite
%   Mixture Models, ch. 6.
%
%   McLachlan's remedy, implemented here: fit both orders to the data and
%   record LR. Then generate B datasets FROM THE FITTED K0 MODEL, refit
%   both orders to each, and record LR_b. The p-value is
%
%       p = (1 + #{LR_b >= LR}) / (B + 1)
%
%   The +1 in both places counts the observed statistic as one draw from
%   the null, which keeps p from ever being exactly 0 (Davison & Hinkley
%   1997, sec. 4.2).
%
%   SIMULATION IS EXACT, not an inversion. In q coordinates each component
%   of the fitted mixture is ALREADY a normalized truncated exponential --
%   that is what removed the normalizer from the likelihood -- so a
%   replicate is drawn by choosing component j with probability q_j and
%   then drawing from that component's truncated law directly: a shifted
%   exponential above xmin in continuous mode, and n_min + Geometric(p_j)
%   with p_j = 1-exp(-lambda_j*dt) in discrete mode. No numerical
%   inversion, no grid, no tail truncation.
%
%   REPLICATES WHERE K1 COLLAPSES ARE KEPT, NOT DISCARDED. Under the null
%   the K1 fit SHOULD usually collapse onto K0, giving LR_b near zero;
%   those replicates are the atom at zero that makes this test different
%   from chi-square, and dropping them would bias the null upward and the
%   p-value down. Only replicates where the fitter found no valid fit at
%   all are dropped, and the count is reported.
%
%   READ NegativeLRFraction BEFORE THE P-VALUE. LR_b cannot be negative in
%   exact arithmetic, since K0 is nested in K1 and the K1 optimum is at
%   least as good. A negative value means the K1 multistart failed to find
%   its own optimum on that replicate. A few percent is normal; a large
%   fraction means the null LRs are systematically too small, which makes
%   p too small -- raise nStartsBase and nStartsPerComponent and re-run.
%
%   NAME-VALUE OPTIONS
%   SamplingInterval  REQUIRED, same units as xmin.
%   DistributionType  "discrete" (default) or "continuous".
%   B                 bootstrap replicates, default 999. A test has a
%                     DECISION BOUNDARY that the estimate must resolve, so
%                     it needs more replicates than an estimate does:
%                     Davison & Hinkley (1997, sec. 4.2) recommend B >= 999
%                     for tests and reserve 100-200 for standard errors.
%                     The Monte Carlo error on the p-value is
%                     sqrt(p(1-p)/B), so a true p of 0.05 carries a 95%
%                     interval of [0.020, 0.080] at B=199 -- the verdict at
%                     alpha=0.05 is then close to a coin flip on noise
%                     alone -- against [0.036, 0.064] at B=999. Both
%                     satisfy the alpha*(B+1)-integer convention, so that
%                     is not what decides it; the sampling error is. Lower
%                     it only for exploratory runs, and say what you used.
%   MaxRate, MinExpectedCount   passed through to the fitter.
%   nStartsBase, nStartsPerComponent, maxStarts
%                     multistart budget for the REPLICATE fits (the fits
%                     to the real data always use the fitter's defaults).
%                     Defaults here are lower than the fitter's, for speed;
%                     raise them if NegativeLRFraction is large.
%   UseParallel       false (default) or true. Runs the replicate loop over
%                     a parallel pool. Each replicate carries its OWN seed,
%                     drawn serially from RandomSeed, so the result is
%                     identical serial or parallel and independent of
%                     worker count -- test 13 asserts exactly that. Needs
%                     no toolbox: without the Parallel Computing Toolbox
%                     MATLAB runs the loop serially anyway. The loop is
%                     embarrassingly parallel with near-equal iteration
%                     cost, so expect close to linear speedup once the pool
%                     has started.
%   RandomSeed, Verbose
%
%   OUTPUT
%   L.K0, L.K1, L.LogLik0, L.LogLik1, L.LR
%   L.pValue          the bootstrap p-value
%   L.B, L.BValid     requested and usable replicates
%   L.LRNull          the simulated null LRs, for plotting or a quantile
%   L.AtomAtZero      fraction of null LRs that are ~0, i.e. how often K1
%                     collapsed onto K0 under the null
%   L.NegativeLRFraction   see above
%   L.pValueChiSquare df and p under the INAPPLICABLE chi-square, reported
%                     only so the size of the error is visible
%   L.Note            plain-language verdict
%   L.n, L.xmin, L.SamplingInterval, L.DistributionType
%
%   EXAMPLE
%       L = HyperexponentialLRT(bouts, 300, 2, 3, SamplingInterval=30);
%       % L.pValue = 0.005  ->  the third component is real
%
%   See also FITHYPEREXPONENTIALMLE, COMPAREBOUTMODELS.

arguments
    eventseries double {mustBeReal}
    xmin (1,1) double {mustBePositive}
    K0 (1,1) double {mustBeInteger,mustBePositive}
    K1 (1,1) double {mustBeInteger,mustBePositive}
    options.SamplingInterval (1,1) double = NaN
    options.DistributionType (1,1) string {mustBeMember(options.DistributionType,["continuous","discrete"])} = "discrete"
    options.B (1,1) double {mustBeInteger,mustBePositive} = 999
    options.MaxRate (1,1) double = NaN
    options.MinExpectedCount (1,1) double {mustBeNonnegative} = 5
    options.nStartsBase (1,1) double {mustBeInteger,mustBePositive} = 4
    options.nStartsPerComponent (1,1) double {mustBeInteger,mustBePositive} = 4
    options.maxStarts (1,1) double {mustBeInteger,mustBePositive} = 20
    options.RandomSeed = []
    options.UseParallel (1,1) logical = false
    options.Verbose (1,1) logical = true
end

if isnan(options.SamplingInterval)
    error('HyperexponentialLRT:SamplingIntervalRequired', ...
        ['SamplingInterval is required, in the SAME UNITS as xmin. It sets ' ...
         'the discrete bin width and n_min, and the replicate simulator ' ...
         'needs it to place draws on the same grid as the data.']);
end
if ~(options.SamplingInterval > 0)
    error('HyperexponentialLRT:InvalidSamplingInterval', ...
        'SamplingInterval must be positive; got %g.', options.SamplingInterval);
end
if K0 >= K1
    error('HyperexponentialLRT:NotNested', ...
        ['K0 must be strictly less than K1 -- the test compares a smaller ' ...
         'mixture against a larger one that CONTAINS it. Got K0=%d, K1=%d.'], ...
        K0, K1);
end
if ~isempty(options.RandomSeed)
    rng(options.RandomSeed);
end

dt = options.SamplingInterval;
isDiscrete = strcmp(options.DistributionType, "discrete");
common = {'SamplingInterval', dt, 'DistributionType', options.DistributionType, ...
          'MinExpectedCount', options.MinExpectedCount, ...
          'ErrorOnNoValidFit', false, 'Verbose', false};
if ~isnan(options.MaxRate), common = [common, {'MaxRate', options.MaxRate}]; end

% ------------------------------------------------- both orders on the data
% MaxComponents=K1 fits 1..K1 in one call, so both orders come from the
% same optimizer settings -- an LR built from two differently-tuned fits
% would confound the comparison with the tuning.
H = FitHyperexponentialMLE(eventseries, xmin, common{:}, 'MaxComponents', K1);
f0 = H.AllFits(K0);
f1 = H.AllFits(K1);
if ~f0.Success
    error('HyperexponentialLRT:NullNotFitted', ...
        ['K0=%d could not be fitted to the data, so there is no null model ' ...
         'to simulate from: %s'], K0, f0.DegenerateReason);
end
if ~f1.Success
    error('HyperexponentialLRT:AlternativeNotFitted', ...
        ['K1=%d could not be fitted to the data: %s'], K1, f1.DegenerateReason);
end
LR = 2 * (f1.LogLik - f0.LogLik);
n = H.n;

if options.Verbose
    fprintf(['HyperexponentialLRT: K=%d vs K=%d on n=%d, %s mode.\n' ...
             '  logL(K=%d) = %.4f,  logL(K=%d) = %.4f,  LR = %.4f\n'], ...
        K0, K1, n, options.DistributionType, K0, f0.LogLik, K1, f1.LogLik, LR);
    if options.UseParallel
        fprintf(['  simulating %d replicates from the fitted K=%d model, ' ...
            'in parallel (per-replicate seeds, so the answer does not ' ...
            'depend on worker count) ...\n'], options.B, K0);
    else
        fprintf('  simulating %d replicates from the fitted K=%d model ...\n', ...
            options.B, K0);
    end
end

% ------------------------------------------------------ the null by simulation
q0 = f0.WeightsObserved(:).';
r0 = f0.Rates(:).';
% Every replicate gets its OWN seed, drawn serially here from the caller's
% stream, and reseeds at the top of its iteration. That makes each replicate
% self-contained, so the null sample no longer depends on the ORDER the
% iterations finish in -- which is what a bare parfor would destroy. The
% same RandomSeed then gives the same p on one core or sixty-four, a
% stronger guarantee than the serial loop had, since that one relied on
% every replicate drawing from one sequential stream.
seeds = randi(2^31-1, 1, options.B);

% parfor's second argument caps the workers; 0 forces the loop to run in
% the client. So one loop body serves both modes, and nothing here needs
% the Parallel Computing Toolbox: without it MATLAB runs parfor as a plain
% loop, and Octave accepts both forms too.
nw = 0;
if options.UseParallel, nw = Inf; end
serial = ~options.UseParallel;
verbose = options.Verbose;
tick = max(1, round(options.B/10));

% Broadcast copies, so the loop body does not reach into `options` (which
% parfor would ship whole to every worker).
B_  = options.B;
nsB = options.nStartsBase;
nsC = options.nStartsPerComponent;
msX = options.maxStarts;

LRb = nan(1, B_);
parfor (b = 1:B_, nw)
    rng(seeds(b));
    try
        sim = simTruncMixture(q0, r0, n, xmin, dt, isDiscrete);
        Hb = FitHyperexponentialMLE(sim, xmin, common{:}, ...
            'MaxComponents', K1, ...
            'nStartsBase', nsB, ...
            'nStartsPerComponent', nsC, ...
            'maxStarts', msX);
        g0 = Hb.AllFits(K0); g1 = Hb.AllFits(K1);
        % Degenerate is NOT a reason to drop a replicate: under the null
        % K1 collapsing onto K0 is the expected outcome and is what puts
        % the atom at zero into the null distribution. A replicate is only
        % lost when no valid fit was found at all, which leaves NaN here
        % and is counted after the loop -- a running counter would be a
        % reduction inside try/catch, which parfor cannot analyse.
        if g0.Success && isfinite(g0.LogLik) && isfinite(g1.LogLik)
            LRb(b) = 2 * (g1.LogLik - g0.LogLik);
        end
    catch
        % leaves NaN
    end
    if serial && verbose && mod(b, tick) == 0
        % Only serially: in a real pool the iterations finish out of order
        % and MATLAB buffers worker output, so per-iteration ticks stop
        % meaning anything.
        fprintf('    %d/%d\n', b, B_);
    end
end
nFailed = sum(~isfinite(LRb));

valid = LRb(isfinite(LRb));
nValid = numel(valid);
if nValid < 20
    error('HyperexponentialLRT:TooFewReplicates', ...
        ['Only %d of %d replicates produced a usable pair of fits, which is ' ...
         'too few for a p-value. Raise the multistart budget or B.'], ...
        nValid, options.B);
end

negFrac = mean(valid < -1e-6);
valid = max(valid, 0);          % see the header: negatives are optimizer noise
LRobs = max(LR, 0);
pBoot = (1 + sum(valid >= LRobs)) / (nValid + 1);
atom  = mean(valid < 1e-6);

% The chi-square that does NOT apply, reported so its error is visible.
dfNaive = 2 * (K1 - K0);
pChi = gammainc(LRobs/2, dfNaive/2, 'upper');

L = struct();
L.K0 = K0; L.K1 = K1;
L.LogLik0 = f0.LogLik; L.LogLik1 = f1.LogLik;
L.LR = LR;
L.pValue = pBoot;
L.B = options.B;
L.UseParallel = options.UseParallel;
L.BValid = nValid;
L.FailedReplicates = nFailed;
L.LRNull = valid;
L.AtomAtZero = atom;
L.NegativeLRFraction = negFrac;
L.pValueChiSquare = pChi;
L.ChiSquareDf = dfNaive;
L.n = n; L.xmin = xmin;
L.SamplingInterval = dt;
L.DistributionType = options.DistributionType;
L.Fit0 = f0; L.Fit1 = f1;

if pBoot <= 0.05
    L.Note = sprintf(['K=%d is supported over K=%d (bootstrap p=%.4g from ' ...
        '%d replicates).'], K1, K0, pBoot, nValid);
else
    L.Note = sprintf(['K=%d is NOT supported over K=%d (bootstrap p=%.4g ' ...
        'from %d replicates); prefer the smaller model.'], K1, K0, pBoot, nValid);
end
if negFrac > 0.1
    L.Note = [L.Note sprintf([' CAUTION: %.0f%% of null LRs came out ' ...
        'negative, so the K=%d replicate fits are under-optimized and this ' ...
        'p-value is too small. Raise nStartsBase and nStartsPerComponent.'], ...
        100*negFrac, K1)];
end

if options.Verbose
    fprintf('  LR = %.4f,  bootstrap p = %.4g  (%d valid replicates)\n', ...
        LRobs, pBoot, nValid);
    fprintf('  null: %.0f%% at zero (K=%d collapsed onto K=%d), 95th pct = %.3f\n', ...
        100*atom, K1, K0, quantileSimple(valid, 0.95));
    fprintf('  for contrast, the INAPPLICABLE chi2(%d) would give p = %.4g\n', ...
        dfNaive, pChi);
    if negFrac > 0
        fprintf('  %.1f%% of null LRs were negative (optimizer noise)\n', 100*negFrac);
    end
    fprintf('  %s\n', L.Note);
end
end

% ========================================================================
function d = simTruncMixture(q, rates, n, xmin, dt, isDiscrete)
%SIMTRUNCMIXTURE Exact draw from a fitted truncated hyperexponential.
% In q coordinates each component is already normalized over the observed
% region, so no rejection and no inversion are needed: pick a component,
% then draw from its own truncated law.
q = q(:).' / sum(q);
rates = rates(:).';
cw = cumsum(q);
u = rand(n, 1);
comp = ones(n, 1);
for j = 1:numel(cw)-1
    comp = comp + (u > cw(j));
end
lam = rates(comp);
lam = lam(:);
v = rand(n, 1);
v(v <= 0) = eps;
if isDiscrete
    n_min = max(1, round(xmin / dt));
    % P(N = n_min + g | component j) = (1-p_j)^g * p_j, p_j = 1-exp(-lam*dt),
    % which is exactly the truncated geometric the discrete likelihood uses.
    log1mp = -lam * dt;                        % log(1-p) without cancellation
    g = floor(log(v) ./ log1mp);
    g(~isfinite(g)) = 0;
    d = (n_min + g) * dt;
else
    % Memorylessness: an exponential conditioned on T >= xmin is xmin plus
    % a fresh exponential of the same rate.
    d = xmin - log(v) ./ lam;
end
end

function q = quantileSimple(x, p)
%QUANTILESIMPLE Order statistic, so no Statistics Toolbox is needed.
x = sort(x(:));
if isempty(x), q = NaN; return; end
idx = max(1, min(numel(x), ceil(p * numel(x))));
q = x(idx);
end
