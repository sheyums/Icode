function R = SleepAmountACF(flies, opts)
%SLEEPAMOUNTACF  Autocorrelation of sleep minutes per window against a renewal null.
%
%   R = SleepAmountACF(flies)
%   R = SleepAmountACF(flies, BinMinutes=30, MaxLagHours=12, ...)
%
% A homeostat predicts COMPENSATION: a window with less sleep than usual is
% followed, hours later, by more. Per fly, the record (from its first bout
% onset to its last bout end) is cut into windows of BinMinutes, the sleep
% minutes in each window counted, and the autocorrelation computed at lags
% 1..MaxLagHours. Binning alone makes neighbouring windows correlate (a long
% bout fills several), so the reference is not zero but a RENEWAL null: the
% same fly's bouts with sleep durations permuted among sleep positions and
% wake among wake positions (RenewalSurrogate), the trace rebuilt and binned
% identically. The deviation curve = observed ACF minus the null mean ACF.
%
% GROUP TEST (max-statistic, corrects for scanning lags): the group-mean
% deviation curve's MINIMUM over lags MinLagHours..MaxLagHours, against the
% same minimum of each surrogate index's group-mean deviation curve;
% one-sided (more negative). Likewise the MAXIMUM (persistence). Per fly:
% its own minimum against its own null minima (one-sided p).
%
% Windows that are not fully covered by bouts (a gap between non-adjacent
% bouts) are NaN and drop out of every pair that uses them.
%
% INPUT   flies: struct array (or cell) with column vectors Dur (s), IsSleep,
%         Onset (s), bouts in time order.
% OPTIONS BinMinutes 30; MaxLagHours 12; MinLagHours 1 (start of the range
%         the extreme is taken over); NumSurrogates 1000; Seed 1 (private
%         'mrg32k3a' stream, substream = fly index -- the SAME surrogates as
%         ProcessSProxy with the same Seed); FlyNames.
% OUTPUT  struct R
%   LagHours        1 x L
%   ACF             nFly x L observed
%   NullMeanACF     nFly x L
%   Deviation       nFly x L
%   GroupDeviation  1 x L (mean over flies)
%   Group           table: MinDev, MinLagHours, MinP, MaxDev, MaxLagHours,
%                   MaxP, NFlies, NFliesMinSig (per-fly min p < 0.05)
%   PerFly          table: Fly, NWindows, MinDev, MinLagHours, MinP
%   GroupNullMin, GroupNullMax   1 x NumSurrogates
%
% MEASURED CALIBRATION (test_SleepAmountACF, 2026-11-06; 48-fly groups of
% SimulateSleepWake flies): independent bouts (100 groups) MinP 0.100 -- ABOVE
% nominal, so the compensation p is not calibrated -- and MaxP 0.030; homeostat
% (50) MinP 1.000 at 1.0 h; shared slow state with no homeostat (50) MinP 0.580.
% Neither calibrated nor specific: use descriptively.
%
% Needs no toolbox.
%
% CHANGE LOG
%   2026-11-06  First version (PRESPEC 2026-11-06 entry C, committed fe307b0).

arguments
    flies
    opts.BinMinutes (1,1) double {mustBePositive} = 30
    opts.MaxLagHours (1,1) double {mustBePositive} = 12
    opts.MinLagHours (1,1) double {mustBeNonnegative} = 1
    opts.NumSurrogates (1,1) double {mustBeInteger, mustBePositive} = 1000
    opts.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.FlyNames string = string.empty
end
if iscell(flies), flies = [flies{:}]; end
nF = numel(flies);
if isempty(opts.FlyNames), names = "fly" + string(1:nF); else, names = opts.FlyNames(:).'; end
bw = opts.BinMinutes * 60;
lags = 1:floor(opts.MaxLagHours * 3600 / bw);
lagH = lags * bw / 3600;
inRange = lagH >= opts.MinLagHours;
nL = numel(lags);
Rn = opts.NumSurrogates;

acf = NaN(nF, nL); nullMean = NaN(nF, nL); nullDev = NaN(nF, nL, Rn); nWin = zeros(nF, 1);
for f = 1:nF
    fl = flies(f);
    acf(f, :) = window_acf(fl, bw, lags);
    nWin(f) = floor((fl.Onset(end) + fl.Dur(end) - fl.Onset(1)) / bw);
    rs = RandStream('mrg32k3a', 'Seed', opts.Seed);
    rs.Substream = f;
    A = NaN(Rn, nL);
    for r = 1:Rn
        A(r, :) = window_acf(RenewalSurrogate(fl, rs), bw, lags);
    end
    nullMean(f, :) = mean(A, 1, 'omitnan');
    nullDev(f, :, :) = permute(A - nullMean(f, :), [3 2 1]);
end
dev = acf - nullMean;

% per fly: own minimum against own null minima
[pfMin, pfIdx] = min(dev(:, inRange), [], 2);
lr = lagH(inRange);
pfP = NaN(nF, 1);
for f = 1:nF
    nm = min(squeeze(nullDev(f, inRange, :)), [], 1);
    pfP(f) = (1 + sum(nm <= pfMin(f))) / (Rn + 1);
end

% group
gDev = mean(dev, 1, 'omitnan');
gNull = squeeze(mean(nullDev, 1, 'omitnan'));        % L x Rn
if nL == 1, gNull = gNull(:).'; end
nullMin = min(gNull(inRange, :), [], 1);
nullMax = max(gNull(inRange, :), [], 1);
[mn, im] = min(gDev(inRange)); [mx, ix] = max(gDev(inRange));
G = table(mn, lr(im), (1 + sum(nullMin <= mn)) / (Rn + 1), mx, lr(ix), ...
    (1 + sum(nullMax >= mx)) / (Rn + 1), nF, sum(pfP < 0.05), ...
    'VariableNames', {'MinDev','MinLagHours','MinP','MaxDev','MaxLagHours','MaxP','NFlies','NFliesMinSig'});

R.LagHours = lagH;
R.ACF = acf; R.NullMeanACF = nullMean; R.Deviation = dev; R.GroupDeviation = gDev;
R.Group = G;
R.PerFly = table(names(:), nWin, pfMin, lr(pfIdx).', pfP, ...
    'VariableNames', {'Fly','NWindows','MinDev','MinLagHours','MinP'});
R.GroupNullMin = nullMin; R.GroupNullMax = nullMax;
R.Settings = opts; R.Settings.FlyNames = names; R.Settings.Version = '2026-11-06';
end

% =====================================================================
function a = window_acf(fl, bw, lags)
y = window_sleep(fl, bw);
a = NaN(1, numel(lags));
for k = 1:numel(lags)
    L = lags(k);
    if L >= numel(y) - 2, continue, end
    u = y(1:end-L); v = y(1+L:end);
    ok = isfinite(u) & isfinite(v);
    if sum(ok) < 3, continue, end
    u = u(ok) - mean(u(ok)); v = v(ok) - mean(v(ok));
    den = sqrt(sum(u.^2) * sum(v.^2));
    if den > 0, a(k) = sum(u .* v) / den; end
end
end

function y = window_sleep(fl, bw)
% Sleep seconds per window, windows from the first onset; NaN where the
% window is not fully covered by bouts.
t0 = fl.Onset(1);
tEnd = fl.Onset(end) + fl.Dur(end);
nW = floor((tEnd - t0) / bw);
e1 = t0 + (0:nW-1).' * bw; e2 = e1 + bw;
a = fl.Onset(:).'; b = a + fl.Dur(:).';
ov = max(0, min(e2, b) - max(e1, a));            % windows x bouts
y = ov * double(fl.IsSleep(:));
covered = sum(ov, 2);
y(covered < bw) = NaN;
end
