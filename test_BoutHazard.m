function test_BoutHazard()
%TEST_BOUTHAZARD Smoke tests for the life-table hazard diagnostic.
%
%   Run from the directory containing BoutHazard.m:
%
%       test_BoutHazard
%
% The point of this function is to decide whether a hazard RISES, and both
% kinds of error matter. Calling a rise on a hyperexponential would excuse
% adding model families that were never needed; missing a real one would
% leave a family excluded that cannot fit. Tests 3 and 4 are those two
% cases, on laws whose hazards are known in closed form.
%
% Sample sizes are large where a band's coverage is being checked, because
% a binomial interval at d/n with n in the hundreds is wide.
%
%  1.  SamplingInterval required     - and non-positive rejected
%  2.  Exponential: hazard is FLAT   - recovered, and the bands cover it
%  3.  Hyperexponential: DECREASING  - no rise called (false-positive guard)
%  4.  Erlang: RISES                 - rise called, bootstrap CI excludes 1
%  5.  The rate conversion           - p/W would be wrong; -log(1-p)/W is not
%  6.  Saturated bins excluded       - the last bin is always d = n
%  7.  Bands bracket the estimate    - and stay non-negative
%  8.  MinAtRisk drops thin bins     - nDropped counts them
%  9.  Truncation is free            - same hazard above a higher xmin
% 10.  RandomSeed reproduces         - identical RiseCI twice

if exist('OCTAVE_VERSION', 'builtin')
    bh = @BoutHazard_oct;
else
    bh = @BoutHazard;
end
nPassed = 0; nFailed = 0;
fprintf('\nRunning tests for BoutHazard...\n\n');

%% Test 1: SamplingInterval is required
try
    d = 100 + (1:500)';
    a = ''; try, bh(d, 100); catch e, a = e.identifier; end
    b = ''; try, bh(d, 100, 'SamplingInterval', 0); catch e, b = e.identifier; end
    ok = ~isempty(strfind(a, 'SamplingIntervalRequired')) && ...
         ~isempty(strfind(b, 'InvalidSamplingInterval')); %#ok<STREMP>
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 1, ...
        'SamplingInterval is required and must be positive', ...
        sprintf('ids "%s" / "%s"', a, b));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 1, '', err.message);
end

%% Test 2: an exponential's hazard is flat at lambda, and the bands cover it
%  The estimator has no idea the law is exponential; it must return lambda
%  across two decades of bins of very unequal width.
try
    rng(1, 'twister');
    lam = 1/500; xm = 100; dt = 1; n = 20000;
    % Engine convention: a bout of k bins lies in ((k-1)dt, k dt], so
    % k = ceil(x/dt) and the shortest attainable duration is xm. round()
    % would put each bout in the NEAREST bin, half a step off, which biases
    % the first bin's rate -- measured at a 7.4% miss rate there against
    % 5% elsewhere.
    d = xm - dt + ceil(-log(rand(n,1))/lam/dt)*dt;
    H = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 12);
    rel = H.Hazard / lam;
    nb = numel(H.Hazard);
    misses = nb - sum(H.Lower <= lam & lam <= H.Upper);
    % A 95% band misses its own parameter 5% of the time BY DESIGN, so the
    % number of misses over nb bins is Binomial(nb, 0.05) -- not zero, and
    % not "at most one". An earlier version asserted at most one miss in
    % 11 bins, which fails for 10.3% of perfectly calibrated datasets
    % (measured over 1000 seeds; P(>=2 misses) = 0.102 for independent
    % bands). MATLAB's seed 1 happened to land there and Octave's did not.
    % The rule is the binomial's own upper tail: fail only above the 99th
    % percentile, so a correct estimator false-fails ~0.15% of the time.
    allowed = binoUpperTail(nb, 0.05, 0.99);
    bad = {};
    if max(abs(rel - 1)) > 0.25
        bad{end+1} = sprintf('hazard off by up to %.0f%%', 100*max(abs(rel-1)));
    end
    if misses > allowed
        bad{end+1} = sprintf(['bands missed lambda in %d of %d bins; ' ...
            'Binomial(%d, 0.05) allows up to %d at the 99th percentile'], ...
            misses, nb, nb, allowed);
    end
    if H.NonMonotone
        bad{end+1} = 'called a rise on a memoryless law';
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 2, ...
        sprintf('flat within %.0f%%, %d of %d bins missed (binomial allows %d)', ...
            100*max(abs(rel-1)), misses, nb, allowed), strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 2, '', err.message);
end

%% Test 3: a hyperexponential's hazard DECREASES -- no rise may be called
%  The false positive that matters. A mixture of exponentials is a sum of
%  decreasing terms, so its hazard decreases at every order; any rise the
%  estimator reports here is binning noise. Binning noise DOES produce an
%  apparent rise -- the assertion is that it is not RESOLVED, not that it
%  is absent.
try
    rng(2, 'twister');
    tau = [60 1500]; w = [0.6 0.4]; xm = 100; dt = 1; n = 20000;
    c = 1 + (rand(n,1) > w(1));
    d = xm + round(-tau(c)' .* log(rand(n,1)));
    H = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 14, ...
           'Bootstrap', 199, 'RandomSeed', 2);
    bad = {};
    if ~(H.Hazard(1) > H.Hazard(end))
        bad{end+1} = 'hazard did not fall across the range';
    end
    if H.NonMonotone
        bad{end+1} = sprintf(['called a rise (x%.2f, CI [%.2f %.2f]) on a ' ...
            'strictly decreasing hazard'], H.RiseRatio, H.RiseCI(1), H.RiseCI(2));
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 3, ...
        sprintf('falls %.3g -> %.3g; apparent rise x%.2f left unresolved', ...
            H.Hazard(1), H.Hazard(end), H.RiseRatio), strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 3, '', err.message);
end

%% Test 4: an Erlang's hazard RISES -- it must be called, with an interval
try
    rng(3, 'twister');
    xm = 100; dt = 1; n = 20000; r = 1/200;
    d = xm + round(-log(rand(n,1))/r - log(rand(n,1))/r - log(rand(n,1))/r);
    H = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 14, ...
           'Bootstrap', 199, 'RandomSeed', 3);
    bad = {};
    if ~H.NonMonotone
        bad{end+1} = sprintf('missed the rise (x%.2f)', H.RiseRatio);
    end
    if ~(H.RiseCI(1) > 1)
        bad{end+1} = sprintf('bootstrap CI [%.2f %.2f] does not exclude 1', ...
            H.RiseCI(1), H.RiseCI(2));
    end
    if H.PeakIndex <= H.TroughIndex
        bad{end+1} = 'peak does not follow the trough';
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 4, ...
        sprintf('rise x%.0f, CI [%.0f %.0f] excludes 1', ...
            H.RiseRatio, H.RiseCI(1), H.RiseCI(2)), strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 4, '', err.message);
end

%% Test 5: the rate conversion, against the arithmetic it replaced
%  -log(1-p)/W recovers lambda exactly; p/W does not, and the error grows
%  with bin width, so on log bins it fakes a DECREASING hazard.
try
    lam = 1/500; xm = 100; dt = 1;
    S = @(t) exp(-lam*(t - xm));
    d = xm + [1; 2; 3];                     % the data only sets the range
    edges = xm - dt + [0, 50, 200, 800, 3000];
    H = bh([d; xm + 3000], xm, 'SamplingInterval', dt, 'BinEdges', edges, ...
           'MinAtRisk', 1, 'SurvivalHandle', S);
    W = H.BinWidth;
    pW = 1 - exp(-lam*W);                   % what p/W would have divided
    naive = pW ./ W;
    bad = {};
    if max(abs(H.ModelHazard - lam)) / lam > 1e-9
        bad{end+1} = sprintf('model hazard off by %.3g', ...
            max(abs(H.ModelHazard - lam))/lam);
    end
    if ~(max(abs(naive - lam))/lam > 0.3)
        bad{end+1} = 'the naive p/W form was not actually wrong here, so this test proves nothing';
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 5, ...
        sprintf('exact to %.1g where p/W errs by %.0f%%', ...
            max(abs(H.ModelHazard-lam))/lam, 100*max(abs(naive-lam))/lam), ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 5, '', err.message);
end

%% Test 6: the saturated final bin is never reported
try
    rng(6, 'twister');
    xm = 10; dt = 1;
    d = xm + round(-50*log(rand(3000,1)));
    H = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 10, 'MinAtRisk', 1);
    bad = {};
    if any(~isfinite(H.Hazard))
        bad{end+1} = 'a non-finite hazard was reported';
    end
    if any(H.Events >= H.AtRisk)
        bad{end+1} = 'a saturated bin (d = n) survived into the output';
    end
    if H.nSaturated < 1
        bad{end+1} = 'no bin was saturated, so this test proves nothing';
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 6, ...
        sprintf('%d saturated bin(s) excluded, all reported rates finite', ...
            H.nSaturated), strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 6, '', err.message);
end

%% Test 7: the bands bracket the estimate and stay non-negative
try
    rng(7, 'twister');
    xm = 100; dt = 1;
    d = xm + round(-300*log(rand(4000,1)));
    H = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 16);
    H9 = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 16, 'Alpha', 0.01);
    bad = {};
    if any(H.Lower > H.Hazard) || any(H.Hazard > H.Upper)
        bad{end+1} = 'estimate outside its own band';
    end
    if any(H.Lower < 0)
        bad{end+1} = 'negative lower band';
    end
    if ~all(H9.Lower <= H.Lower + 1e-12) || ~all(H9.Upper >= H.Upper - 1e-12)
        bad{end+1} = 'the 99% band is not wider than the 95% band';
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 7, ...
        'bands bracket the estimate, stay >= 0, and widen with 1-Alpha', ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 7, '', err.message);
end

%% Test 8: MinAtRisk drops thin bins and counts them
try
    rng(8, 'twister');
    xm = 100; dt = 1;
    d = xm + round(-300*log(rand(2000,1)));
    Hlo = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 24, 'MinAtRisk', 1);
    Hhi = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 24, 'MinAtRisk', 200);
    bad = {};
    if ~(numel(Hhi.Hazard) < numel(Hlo.Hazard))
        bad{end+1} = 'a higher MinAtRisk did not drop any bin';
    end
    if any(Hhi.AtRisk < 200)
        bad{end+1} = 'a bin below MinAtRisk survived';
    end
    if Hhi.nDropped ~= (numel(Hhi.BinEdges) - 1 - numel(Hhi.Hazard))
        bad{end+1} = 'nDropped does not account for the missing bins';
    end
    [nPassed, nFailed] = rep(isempty(bad), nPassed, nFailed, 8, ...
        sprintf('%d bins at MinAtRisk=1, %d at 200, nDropped=%d', ...
            numel(Hlo.Hazard), numel(Hhi.Hazard), Hhi.nDropped), ...
        strjoin(bad, '; '));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 8, '', err.message);
end

%% Test 9: truncation is free -- the hazard above a higher xmin is the same
%  The claim that makes this the right statistic for left-truncated bouts.
%  Estimating from bouts >= 300 must agree with estimating from bouts
%  >= 100 and reading the same region, because the hazard is conditional.
try
    rng(9, 'twister');
    dt = 1; n = 40000;
    d = 100 + round(-400*log(rand(n,1)));
    edges = 299 + [0, 200, 500, 1200, 2500];
    Ha = bh(d, 100, 'SamplingInterval', dt, 'BinEdges', edges, 'MinAtRisk', 20);
    Hb = bh(d(d >= 300), 300, 'SamplingInterval', dt, 'BinEdges', edges, 'MinAtRisk', 20);
    rel = abs(Ha.Hazard - Hb.Hazard) ./ Ha.Hazard;
    ok = numel(Ha.Hazard) == numel(Hb.Hazard) && max(rel) < 1e-12;
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 9, ...
        sprintf('identical over %d shared bins (max rel %.1g)', ...
            numel(Ha.Hazard), max(rel)), ...
        sprintf('differs by up to %.3g', max(rel)));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 9, '', err.message);
end

%% Test 10: RandomSeed reproduces the bootstrap interval
try
    rng(10, 'twister');
    xm = 100; dt = 1;
    d = xm + round(-log(rand(3000,1))/0.004 - log(rand(3000,1))/0.004);
    A = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 12, 'Bootstrap', 99, 'RandomSeed', 42);
    B = bh(d, xm, 'SamplingInterval', dt, 'NumBins', 12, 'Bootstrap', 99, 'RandomSeed', 42);
    ok = isequal(A.RiseCI, B.RiseCI) && all(isfinite(A.RiseCI));
    [nPassed, nFailed] = rep(ok, nPassed, nFailed, 10, ...
        sprintf('same seed, same CI [%.3f %.3f]', A.RiseCI(1), A.RiseCI(2)), ...
        sprintf('[%.4f %.4f] vs [%.4f %.4f]', A.RiseCI, B.RiseCI));
catch err
    [nPassed, nFailed] = rep(false, nPassed, nFailed, 10, '', err.message);
end

fprintf('\nSummary: %d passed, %d failed, %d total\n', ...
    nPassed, nFailed, nPassed + nFailed);
end

function k = binoUpperTail(nb, p, conf)
%BINOUPPERTAIL Smallest k with P(X <= k) >= conf for X ~ Binomial(nb, p).
% Computed from the pmf rather than binoinv, which needs the Statistics
% Toolbox -- this suite runs without one.
c = 0; k = nb;
for i = 0:nb
    c = c + nchoosek(nb, i) * p^i * (1-p)^(nb-i);
    if c >= conf, k = i; return; end
end
end

function [nP, nF] = rep(ok, nP, nF, idx, passMsg, failMsg)
if ok
    fprintf('[PASS] Test %d: %s\n', idx, passMsg); nP = nP + 1;
else
    fprintf('[FAIL] Test %d: %s\n', idx, failMsg); nF = nF + 1;
end
end
