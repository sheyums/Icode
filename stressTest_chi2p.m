function stressTest_chi2p
%STRESSTEST_CHI2P Stress tests for chi2p.
%
%   STRESSTEST_CHI2P exercises the chi-square periodogram and, more
%   importantly, its significance machinery. A periodogram that reports the
%   right period but the wrong p-value is worse than useless in a screen,
%   so most of these tests are about calibration rather than about the
%   statistic itself.
%
%   The suite evaluates:
%
%       1. Input validation and option plumbing
%       2. Period recovery for known synthetic rhythms
%       3. Waveform independence: square, sinusoid, and spiked profiles
%       4. Null calibration under the block-permutation model
%       5. The chi-square null is anti-conservative on autocorrelated data
%       6. Power against strong and moderate rhythms
%       7. Reproducibility and determinism
%       8. Degenerate and boundary inputs
%
%   GROUP 5 IS A REGRESSION GUARD, NOT A BUG
%   ----------------------------------------
%   It asserts that NullModel="chi2" over-rejects on autocorrelated input.
%   That is a property of the textbook test, and the reason the default is
%   the permutation null. The test exists so that nobody "simplifies"
%   chi2p by making the chi-square threshold the default again. Group 4
%   pins the behaviour that matters: the default null must hold its
%   nominal error rate on the SAME autocorrelated input.
%
%   Synthetic autocorrelated records are used rather than real data so the
%   suite is self-contained. They are built by convolving Poisson counts
%   with a short kernel, which reproduces the bout structure responsible
%   for breaking the chi-square assumption.
%
%   COMPATIBILITY
%   -------------
%   Requires chi2p.m on the MATLAB path.
%
%   See also CHI2P, FITPERIODICLOCOMOTORMODEL

clc;
if isempty(which('chi2p'))
    error('stressTest_chi2p:FunctionNotFound', 'chi2p.m must be on the path.');
end

passed = 0; failed = 0; total = 0; failures = strings(0,1);

% Keep the suite quick: short records, few surrogates. Power estimates are
% therefore coarse, and the thresholds below are set accordingly loose.
NS      = 60;
BIN     = minutes(30);
STEP    = seconds(1);
BAND    = hours([16 32]);
DAYS    = 8;
% The generator below produces ~3 h bouts, so the default 2 h BlockLength
% would be mis-specified for these records and every call would (correctly)
% raise chi2p:BlockShorterThanAutocorrelation. Use an adequate block
% throughout; Group 4E deliberately uses a short one to exercise that path.
BLOCK   = hours(12);

fprintf('====================================================\n');
fprintf('chi2p - stress test\n');
fprintf('====================================================\n\n');

%% Group 1: Input validation
fprintf('-- Group 1: Input validation ------------------------\n');
runErrorTest("1A: empty activity",        @() chi2p([]));
runErrorTest("1B: NaN in activity",       @() chi2p([1 2 NaN 4]));
runErrorTest("1C: Inf in activity",       @() chi2p([1 2 Inf 4]));
runErrorTest("1D: complex activity",      @() chi2p([1 2 3i 4]));
runErrorTest("1E: matrix activity",       @() chi2p(ones(4,4)));
runErrorTest("1F: reversed PeriodRange",  @() chi2p(ones(1000,1), PeriodRange=hours([32 16])));
runErrorTest("1G: bins finer than samples", ...
    @() chi2p(ones(1000,1), InputInterval=minutes(10), BinInterval=minutes(1)));
runErrorTest("1H: period below Nyquist", ...
    @() chi2p(ones(1e5,1), BinInterval=hours(10), PeriodRange=hours([16 32])));
runErrorTest("1I: Alpha >= 1",            @() chi2p(ones(1000,1), Alpha=1));
runErrorTest("1J: zero surrogates",       @() chi2p(ones(1000,1), NumSurrogates=0));
runErrorTest("1K: bad NullModel",         @() chi2p(ones(1000,1), NullModel="bootstrap"));
fprintf('\n');

%% Group 2: Period recovery
fprintf('-- Group 2: Period recovery -------------------------\n');
for tauTrue = [20 24 24.5 28]
    x = makeGated(tauTrue, DAYS, STEP, 0);
    r = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, ...
        NumSurrogates=NS, RandomSeed=1);
    err = r.TauHours - tauTrue;
    runConditionTest(sprintf("2: tau=%g h recovered within one bin", tauTrue), ...
        r.IsRhythmic && abs(err) <= hours(BIN) + 1e-9, ...
        sprintf('got %.2f h (error %+.2f h), p=%.4f', r.TauHours, err, r.PValue));
end
fprintf('\n');

%% Group 3: Waveform independence
fprintf('-- Group 3: Waveform independence -------------------\n');
tau = 24;
shapes = struct( ...
    'name', {"square (12h on/off)", "sinusoid", "narrow spike (2h/cycle)"}, ...
    'fn',   {@(t) double(mod(t,tau) < tau/2), ...
             @(t) 1 + cos(2*pi*t/tau), ...
             @(t) double(mod(t,tau) < 2)});
for s = 1:numel(shapes)
    t = (0:seconds(STEP):DAYS*24*3600-1).'/3600;
    x = round(100*shapes(s).fn(t));
    r = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, ...
        NumSurrogates=NS, RandomSeed=1);
    runConditionTest("3: " + shapes(s).name + " detected at tau=24 h", ...
        r.IsRhythmic && abs(r.TauHours - tau) <= hours(BIN) + 1e-9, ...
        sprintf('tau=%.2f p=%.4f', r.TauHours, r.PValue));
end
fprintf('\n');

%% Group 4: Null calibration, block-permutation model
fprintf('-- Group 4: Null calibration (the critical group) ---\n');
% BlockLength must exceed the autocorrelation time. The generator produces
% 3 h bouts, so 12 h blocks are used here; see Group 4E for what happens
% when that constraint is violated.
warnState = warning('off','chi2p:BlockShorterThanAutocorrelation');
nRec = 40; rejects = 0;
for k = 1:nRec
    x = makeAutocorrNoise(DAYS, STEP, k);
    r = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, ...
        BlockLength=BLOCK, NumSurrogates=NS, RandomSeed=k, Alpha=0.05);
    rejects = rejects + r.IsRhythmic;
end
rate = rejects/nRec;
% With 40 records and NS=60 the Monte Carlo error is large; require only
% that the rate is in the same neighbourhood as nominal, not exactly 0.05.
runConditionTest("4A: permutation null holds near nominal 5% with adequate blocks", ...
    rate <= 0.15, sprintf('rejection rate %.3f over %d records', rate, nRec));
runConditionTest("4B: p-values are not degenerate", ...
    rate < 1.0, sprintf('rejection rate %.3f', rate));

x = makeAutocorrNoise(DAYS, STEP, 99);
r = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, NumSurrogates=NS);
runConditionTest("4C: p-value respects the surrogate resolution floor", ...
    r.PValue >= 1/(NS+1) - 1e-12, sprintf('p=%.5f floor=%.5f', r.PValue, 1/(NS+1)));
runConditionTest("4D: p-value never exceeds one", r.PValue <= 1, sprintf('p=%.4f', r.PValue));

% 4E documents the failure mode that motivates the guard: blocks shorter
% than the autocorrelation inflate the rejection rate.
nRec = 40; rejShort = 0;
for k = 1:nRec
    x = makeAutocorrNoise(DAYS, STEP, k);
    rr = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, ...
        BlockLength=hours(1), NumSurrogates=NS, RandomSeed=k);
    rejShort = rejShort + rr.IsRhythmic;
end
runConditionTest("4E: too-short blocks are anti-conservative (documented failure)", ...
    rejShort/nRec > rate, sprintf('1 h blocks %.3f vs 12 h blocks %.3f', ...
    rejShort/nRec, rate));
warning(warnState);

% 4F: the guard must actually fire for the case 4E just demonstrated.
lastwarn('','');
chi2p(makeAutocorrNoise(DAYS, STEP, 1), InputInterval=STEP, BinInterval=BIN, ...
    PeriodRange=BAND, BlockLength=hours(1), NumSurrogates=5);
[~, wid] = lastwarn;
runConditionTest("4F: guard warns when BlockLength is too short", ...
    strcmp(wid,'chi2p:BlockShorterThanAutocorrelation'), sprintf('got "%s"', wid));

lastwarn('','');
warning('off','chi2p:BlockShorterThanAutocorrelation');
rAc = chi2p(makeAutocorrNoise(DAYS, STEP, 1), InputInterval=STEP, BinInterval=BIN, ...
    PeriodRange=BAND, BlockLength=hours(12), NumSurrogates=5);
[~, wid2] = lastwarn;
warning(warnState);
runConditionTest("4G: guard stays quiet when BlockLength is adequate", ...
    ~strcmp(wid2,'chi2p:BlockShorterThanAutocorrelation'), sprintf('got "%s"', wid2));
runConditionTest("4H: autocorrelation time is reported and positive", ...
    isfinite(rAc.AutocorrelationTimeHours) && rAc.AutocorrelationTimeHours > 0, ...
    sprintf('%g h', rAc.AutocorrelationTimeHours));
fprintf('     (estimated autocorrelation time %.2f h for 3 h synthetic bouts)\n', ...
    rAc.AutocorrelationTimeHours);
fprintf('\n');

%% Group 5: The chi-square null over-rejects (regression guard)
fprintf('-- Group 5: chi2 null is anti-conservative ----------\n');
nRec = 30; rejChi = 0; rejPerm = 0;
for k = 1:nRec
    x = makeAutocorrNoise(DAYS, STEP, 200+k);
    rc = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, NullModel="chi2");
    rp = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, ...
        NullModel="blockPermutation", NumSurrogates=NS, RandomSeed=k);
    rejChi  = rejChi  + rc.IsRhythmic;
    rejPerm = rejPerm + rp.IsRhythmic;
end
% Both bounds matter. The first says the chi-square null is badly
% anti-conservative on this input; the second says the permutation null is
% not, on the SAME input. If the generator's bout length is ever shortened
% below BinInterval the first assertion will fail, which is the intended
% warning that the test has stopped exercising the regime it was written for.
runConditionTest("5A: chi2 null over-rejects on autocorrelated input", ...
    rejChi >= 0.15*nRec, sprintf('chi2 rejected only %d/%d (%.0f%%)', ...
    rejChi, nRec, 100*rejChi/nRec));
runConditionTest("5A2: permutation null does not, on identical records", ...
    rejPerm <= 0.20*nRec && rejChi > rejPerm, ...
    sprintf('chi2 %d/%d vs permutation %d/%d', rejChi, nRec, rejPerm, nRec));
fprintf('     (chi2 %d/%d, permutation %d/%d -- see the header note)\n', ...
    rejChi, nRec, rejPerm, nRec);
% The two nulls differ only in how significance is judged; the underlying
% periodogram, and therefore the peak, must be identical.
xc = makeGated(24, DAYS, STEP, 11);
rc = chi2p(xc, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, NullModel="chi2");
rp = chi2p(xc, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, ...
    NullModel="blockPermutation", NumSurrogates=NS, RandomSeed=1);
runConditionTest("5B: both nulls report the same peak period", ...
    isequaln(rc.TauHours, rp.TauHours), ...
    sprintf('chi2 %.2f h vs permutation %.2f h', rc.TauHours, rp.TauHours));
runConditionTest("5C: both nulls report the same periodogram", ...
    isequaln(rc.Periodogram.Q, rp.Periodogram.Q), 'periodograms differ');
fprintf('\n');

%% Group 6: Power
fprintf('-- Group 6: Power against real-shaped rhythms -------\n');
nRec = 20; det = 0;
for k = 1:nRec
    x = makeAutocorrNoise(DAYS, STEP, 300+k);
    t = (0:numel(x)-1).'*seconds(STEP)/3600;
    x = x .* double(mod(t,24.2) < 12.1);        % strong consolidated rhythm
    r = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, ...
        NumSurrogates=NS, RandomSeed=k);
    det = det + r.IsRhythmic;
end
runConditionTest("6A: >=70% power against a gated 12h on/off rhythm", ...
    det/nRec >= 0.70, sprintf('detected %d/%d', det, nRec));

nRec = 20; det = 0;
for k = 1:nRec
    x = makeAutocorrNoise(DAYS, STEP, 400+k);
    r = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, ...
        NumSurrogates=NS, RandomSeed=k);
    det = det + r.IsRhythmic;
end
runConditionTest("6B: same records without a rhythm are mostly not detected", ...
    det/nRec <= 0.25, sprintf('detected %d/%d', det, nRec));
fprintf('\n');

%% Group 7: Reproducibility
fprintf('-- Group 7: Reproducibility -------------------------\n');
x = makeGated(24, DAYS, STEP, 7);
r1 = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, ...
    NumSurrogates=NS, RandomSeed=42);
r2 = chi2p(x, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, ...
    NumSurrogates=NS, RandomSeed=42);
runConditionTest("7A: identical seed gives identical p-value", ...
    isequaln(r1.PValue, r2.PValue), sprintf('%.6f vs %.6f', r1.PValue, r2.PValue));
runConditionTest("7B: identical seed gives identical tau", ...
    isequaln(r1.TauHours, r2.TauHours), '');
runConditionTest("7C: does not disturb the global RNG stream", ...
    localRngUntouched(x, STEP, BIN, BAND, NS), 'global rng state changed');
runConditionTest("7D: periodogram is seed-independent", ...
    isequaln(r1.Periodogram.Q, r2.Periodogram.Q), '');
fprintf('\n');

%% Group 8: Degenerate and boundary inputs
fprintf('-- Group 8: Degenerate and boundary inputs ----------\n');
xConst = ones(DAYS*24*3600,1);
r = chi2p(xConst, InputInterval=STEP, BinInterval=BIN, PeriodRange=BAND, BlockLength=BLOCK, NumSurrogates=NS);
runConditionTest("8A: constant record is not rhythmic", ~r.IsRhythmic, '');
runConditionTest("8B: constant record returns NaN tau", isnan(r.TauHours), ...
    sprintf('%g', r.TauHours));
runConditionTest("8C: constant record returns NaN p", isnan(r.PValue), '');

% 8D and 8E use deliberately short records to probe structural edges, where
% no BlockLength can be both adequate and leave enough blocks to permute.
% The guard correctly objects; it is not what these two tests are about.
warnState8 = warning('off','chi2p:BlockShorterThanAutocorrelation');

% Two complete cycles is the minimum evidence the statistic will accept.
xShort = makeGated(24, 2, STEP, 0);
r = chi2p(xShort, InputInterval=STEP, BinInterval=BIN, PeriodRange=hours([20 26]), ...
    BlockLength=hours(6), NumSurrogates=NS, RandomSeed=1);
runConditionTest("8D: two-cycle record still produces a finite statistic", ...
    ~isnan(r.Ratio), sprintf('ratio=%g', r.Ratio));

% Periods needing more than the available record must be skipped, not error.
r = chi2p(makeGated(24, 3, STEP, 0), InputInterval=STEP, BinInterval=BIN, ...
    PeriodRange=hours([16 40]), NumSurrogates=NS, RandomSeed=1);
longOnes = r.Periodogram.TauHours > 36;
runConditionTest("8E: periods without two cycles are NaN, not errors", ...
    all(isnan(r.Periodogram.Q(longOnes))), '');
warning(warnState8);
runConditionTest("8F: row-vector input accepted", ...
    ~isempty(chi2p(makeGated(24,4,STEP,0).', InputInterval=STEP, BinInterval=BIN, ...
        PeriodRange=BAND, BlockLength=BLOCK, NumSurrogates=5)), '');
runErrorTest("8G: BlockLength as long as the record", ...
    @() chi2p(makeGated(24,3,STEP,0), InputInterval=STEP, BinInterval=BIN, ...
        PeriodRange=BAND, BlockLength=hours(100), NumSurrogates=5));
runConditionTest("8H: BinMethod mean runs and agrees on tau", ...
    localBinMethodAgrees(STEP, BIN, BAND, NS), 'sum and mean disagree on tau');
fprintf('\n');

%% Summary
fprintf('====================================================\n');
fprintf('SUMMARY: %d passed | %d failed | %d total\n', passed, failed, total);
fprintf('====================================================\n');
if failed > 0
    fprintf('\nFailed tests:\n');
    for k = 1:numel(failures)
        fprintf('  %s\n', failures(k));
    end
    error('stressTest_chi2p:TestsFailed', '%d stress test(s) failed.', failed);
else
    fprintf('\nAll stress tests passed.\n');
end

%% ---------------- nested helpers (share the counters) ----------------
    function runConditionTest(name, cond, note)
        total = total + 1;
        if isscalar(cond) && islogical(cond) && cond
            passed = passed + 1;
            fprintf('  [PASS] %s\n', name);
        else
            failed = failed + 1;
            failures(end+1,1) = name + " -- " + string(note);
            fprintf('  [FAIL] %s -- %s\n', name, note);
        end
    end

    function runErrorTest(name, fh)
        total = total + 1;
        try
            fh();
            failed = failed + 1;
            failures(end+1,1) = name + " -- did not error";
            fprintf('  [FAIL] %s -- did not error\n', name);
        catch
            passed = passed + 1;
            fprintf('  [PASS] %s\n', name);
        end
    end
end

%% ---------------- local generators ----------------
function x = makeGated(tauHours, days, step, seed)
%MAKEGATED Autocorrelated counts confined to the first half of each cycle.
x = makeAutocorrNoise(days, step, seed);
t = (0:numel(x)-1).'*seconds(step)/3600;
x = x .* double(mod(t, tauHours) < tauHours/2);
end

function x = makeAutocorrNoise(days, step, seed)
%MAKEAUTOCORRNOISE Poisson counts smeared into bouts.
%   Plain Poisson noise is independent and would not exercise the very
%   assumption these tests exist to check, so the counts are convolved with
%   a short kernel to create bout-scale autocorrelation.
% Drawn from an explicit RandStream so the suite is reproducible and never
% disturbs the caller's global RNG. Bernoulli rather than Poisson keeps this
% free of any toolbox dependency.
%
% BOUT LENGTH MATTERS. The autocorrelation that breaks the chi-square
% assumption must survive binning: bouts shorter than BinInterval are
% averaged away and the record behaves as if independent, at which point
% the chi-square null is perfectly well calibrated and Group 5 has nothing
% to detect. Measured rejection rates for the chi-square null against
% 30-minute bins: 2 min bouts 7%, 30 min 3%, 90 min 20%, 180 min 33%. The
% 180-minute value reproduces what real per0 records show (about 35%), so
% that is what is used here. Consolidated multi-hour rest and activity
% episodes are ordinary in Drosophila, so this is not a contrived regime.
rs = RandStream('twister','Seed', 1000+seed);
n = round(days*24*3600/seconds(step));
raw = double(rand(rs, n, 1) < 0.02);
kern = ones(round(180*60/seconds(step)),1);      % ~3 h bouts
x = conv(raw, kern, 'same');
end

function tf = localRngUntouched(x, step, bin, band, ns)
rng(12345);
before = rng;
chi2p(x, InputInterval=step, BinInterval=bin, PeriodRange=band, ...
    BlockLength=hours(12), NumSurrogates=ns);
after = rng;
tf = isequal(before.State, after.State) && isequal(before.Type, after.Type);
end

function tf = localBinMethodAgrees(step, bin, band, ns)
x = makeGated(24, 8, step, 3);
rs = chi2p(x, InputInterval=step, BinInterval=bin, PeriodRange=band, ...
    BlockLength=hours(12), NumSurrogates=ns, RandomSeed=1, BinMethod="sum");
rm = chi2p(x, InputInterval=step, BinInterval=bin, PeriodRange=band, ...
    BlockLength=hours(12), NumSurrogates=ns, RandomSeed=1, BinMethod="mean");
tf = isequaln(rs.TauHours, rm.TauHours);
end
