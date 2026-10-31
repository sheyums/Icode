function R = rise_test_20r(x, xmin, label, nullReps)
% RISE_TEST_20R  The 2026-09-15 rise test, unchanged: null = FitHyperexponentialMLE
% MaxComponents 2, discrete, SamplingInterval 1, RandomSeed 1, its K = 2
% SurvivalHandle; BoutHazard at NumBins {12,16,20,24}, Spacing "log",
% MinAtRisk 10, Alpha 0.05, Bootstrap 0, NullReplicates 999, RandomSeed 1.
% Verdict by the 2026-09-15 rule. Run BY NAME from the Icode folder.
%
%   x     : bout durations in seconds (NaN ignored)
%   xmin  : truncation point (2 s wake, 300 s sleep)
%   label : text for the printout
%   nullReps : NullReplicates (default 999; the power follow-up uses 199)

assert(isequal(fileparts(which('BoutHazard')), pwd), ...
    'Run from the Icode folder: BoutHazard resolves to %s', which('BoutHazard'));
assert(isequal(fileparts(which('FitHyperexponentialMLE')), pwd), ...
    'FitHyperexponentialMLE resolves to %s', which('FitHyperexponentialMLE'));

if nargin < 4 || isempty(nullReps), nullReps = 999; end
x = x(isfinite(x) & x >= xmin); x = x(:);
F  = FitHyperexponentialMLE(x, xmin, MaxComponents=2, DistributionType="discrete", ...
    SamplingInterval=1, RandomSeed=1, Verbose=false);
f2 = F.AllFits(2);

nb = [12 16 20 24];
R = struct('label', label, 'n', numel(x), 'tau', f2.Tau(:)', 'q', f2.WeightsObserved(:)', ...
    'logL', f2.LogLik, 'Degenerate', f2.Degenerate, 'NumBins', nb);
for i = 1:4
    H = BoutHazard(x, xmin, SamplingInterval=1, NumBins=nb(i), Spacing="log", MinAtRisk=10, ...
        Alpha=0.05, Bootstrap=0, NullSurvivalHandle=f2.SurvivalHandle, NullReplicates=nullReps, ...
        RandomSeed=1);
    R.trough(i)    = H.BinCenters(H.TroughIndex);
    R.peak(i)      = H.BinCenters(H.PeakIndex);
    R.RiseRatio(i) = H.RiseRatio;
    R.nullMedian(i) = H.RiseNullQuantiles(1);
    R.null95(i)    = H.RiseNullQuantiles(2);
    R.RiseNullP(i) = H.RiseNullP;
    R.NullDropped(i) = H.NullDropped;
    R.nDropped(i)  = H.nDropped; R.nSaturated(i) = H.nSaturated;
end

if any(R.NullDropped > 0)
    R.verdict = 'not computable (NullDropped > 0)';
elseif all(R.RiseNullP < 0.05)
    R.verdict = 'supported';
elseif any(R.RiseNullP < 0.05)
    R.verdict = 'supported but bin-sensitive';
else
    R.verdict = 'not supported';
end

fprintf('%s: n = %d, xmin %g s\n', label, R.n, xmin);
fprintf('  null K=2: tau %.2f %.2f s, q %.4f %.4f, logL %.4f, degenerate %d\n', ...
    R.tau, R.q, R.logL, R.Degenerate);
for i = 1:4
    fprintf(['  bins %2d: trough %7.1f s, peak %7.1f s, RiseRatio %.3f, null median %.3f, ' ...
        '95th %.3f, RiseNullP %.3f, NullDropped %d (bins dropped %d, saturated %d)\n'], ...
        nb(i), R.trough(i), R.peak(i), R.RiseRatio(i), R.nullMedian(i), R.null95(i), ...
        R.RiseNullP(i), R.NullDropped(i), R.nDropped(i), R.nSaturated(i));
end
fprintf('  verdict: %s\n', R.verdict);
end
