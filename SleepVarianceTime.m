function R = SleepVarianceTime(flies, opts)
%SLEEPVARIANCETIME  Is cumulative sleep mean-reverting? Variance-time against a renewal null.
%
%   R = SleepVarianceTime(flies)
%   R = SleepVarianceTime(flies, ScalesHours=[1 2 4 8 16 24], ...)
%
% A homeostat regulates HOW MUCH sleep accumulates. Let Y(t) be the fly's
% cumulative sleep surplus: time asleep since the record start minus pbar*t,
% pbar its overall sleep fraction. Over a stretch of L hours the surplus
% changes by D_L(t) = (sleep in [t, t+L)) - pbar*L. Without regulation (the
% same bouts in random order: a renewal process) Var(D_L) grows about in
% proportion to L at long L, like a random walk. If deficits are repaid,
% Var(D_L) grows more slowly -- the curve bends DOWN (mean reversion). A
% slowly varying state that lengthens bouts for a while bends it UP
% (persistence).
%
% Per fly and scale L: V(L) = variance of D_L(t) over start times t on a
% StepMinutes grid (overlapping windows), and the ratio to the renewal null
%   logR(L) = log( V(L) / mean over surrogates of V_surrogate(L) ).
% Per fly, Slope = least-squares slope of logR(L) on log2(L) over the
% scales: negative = mean reversion growing with timescale, positive =
% persistence. GROUP: mean Slope over flies, against the same mean on each
% surrogate index (surrogate r's own logR uses the same null mean);
% one-sided p for reversion (more negative) and for persistence.
%
% RUNS (bout level, secondary): per fly and bout type, the number of runs of
% bouts above/below the type's median (bouts equal to the median dropped);
% more runs than the renewal surrogates = reversion, fewer = persistence.
% Z = (observed - surrogate mean)/surrogate SD; group mean Z against the
% surrogates' group means, two-sided.
%
% The surrogates are RenewalSurrogate draws: with the same Seed, identical to
% those of SleepAmountACF and ProcessSProxy.
%
% INPUT   flies: struct array (or cell) with Dur (s), IsSleep, Onset (s), in
%         time order, CONTIGUOUS (each bout starts where the previous ends;
%         an error otherwise).
% OPTIONS ScalesHours [1 2 4 8 16 24]; StepMinutes 15; NumSurrogates 1000;
%         Seed 1; FlyNames.
% OUTPUT  struct R
%   ScalesHours
%   LogR        nFly x nScale
%   Slope       nFly x 1
%   Group       table: MeanSlope, P_Reversion, P_Persistence, MeanLogR (one
%               column per scale as a vector), RunsZ_Sleep, RunsP_Sleep,
%               RunsZ_Wake, RunsP_Wake, NFlies
%   PerFly      table: Fly, RecordHours, Slope, Slope_P_Reversion,
%               RunsZ_Sleep, RunsZ_Wake
%
% MEASURED CALIBRATION: see test_SleepVarianceTime.
%
% Needs no toolbox.
%
% CHANGE LOG
%   2026-11-06  First version (PRESPEC 2026-11-06 entry E).

arguments
    flies
    opts.ScalesHours (1,:) double {mustBePositive} = [1 2 4 8 16 24]
    opts.StepMinutes (1,1) double {mustBePositive} = 15
    opts.NumSurrogates (1,1) double {mustBeInteger, mustBePositive} = 1000
    opts.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.FlyNames string = string.empty
end
if iscell(flies), flies = [flies{:}]; end
nF = numel(flies);
if isempty(opts.FlyNames), names = "fly" + string(1:nF); else, names = opts.FlyNames(:).'; end
L = opts.ScalesHours * 3600; nL = numel(L);
x = log2(opts.ScalesHours(:)); xc = x - mean(x);
Rn = opts.NumSurrogates;

logR = NaN(nF, nL); slope = NaN(nF, 1); nullSlope = NaN(nF, Rn); pf = NaN(nF, 1);
rz = NaN(nF, 2); rzNull = NaN(nF, 2, Rn); recH = NaN(nF, 1);
for f = 1:nF
    fl = flies(f);
    if any(fl.Onset(2:end) ~= fl.Onset(1:end-1) + fl.Dur(1:end-1))
        error('SleepVarianceTime:Gap', 'Fly %d is not contiguous.', f);
    end
    recH(f) = (fl.Onset(end) + fl.Dur(end) - fl.Onset(1)) / 3600;
    Vo = vt(fl, L, opts.StepMinutes * 60);
    ro = runs(fl);
    rs = RandStream('mrg32k3a', 'Seed', opts.Seed);
    rs.Substream = f;
    Vn = NaN(Rn, nL); rn = NaN(Rn, 2);
    for r = 1:Rn
        s = RenewalSurrogate(fl, rs);
        Vn(r, :) = vt(s, L, opts.StepMinutes * 60);
        rn(r, :) = runs(s);
    end
    mV = mean(Vn, 1, 'omitnan');
    logR(f, :) = log(Vo ./ mV);
    ok = isfinite(logR(f, :));
    if sum(ok) >= 3
        slope(f) = sum(xc(ok).' .* (logR(f, ok) - mean(logR(f, ok)))) / sum(xc(ok).^2);
        lr = log(Vn ./ mV);
        for r = 1:Rn
            y = lr(r, ok); nullSlope(f, r) = sum(xc(ok).' .* (y - mean(y))) / sum(xc(ok).^2);
        end
        pf(f) = (1 + sum(nullSlope(f, :) <= slope(f))) / (Rn + 1);
    end
    sd = std(rn, 0, 1, 'omitnan');
    rz(f, :) = (ro - mean(rn, 1, 'omitnan')) ./ sd;
    rzNull(f, :, :) = permute((rn - mean(rn, 1, 'omitnan')) ./ sd, [3 2 1]);
end
in = isfinite(slope);
ms = mean(slope(in)); nd = mean(nullSlope(in, :), 1, 'omitnan');
zs = mean(rz(in, :), 1, 'omitnan'); zn = squeeze(mean(rzNull(in, :, :), 1, 'omitnan'));
twoP = @(o, n) min(1, 2 * min(1 + sum(n >= o), 1 + sum(n <= o)) / (numel(n) + 1));
G = table(ms, (1 + sum(nd <= ms)) / (Rn + 1), (1 + sum(nd >= ms)) / (Rn + 1), ...
    mean(logR(in, :), 1), zs(1), twoP(zs(1), zn(1, :)), zs(2), twoP(zs(2), zn(2, :)), sum(in), ...
    'VariableNames', {'MeanSlope','P_Reversion','P_Persistence','MeanLogR', ...
    'RunsZ_Sleep','RunsP_Sleep','RunsZ_Wake','RunsP_Wake','NFlies'});
R.ScalesHours = opts.ScalesHours; R.LogR = logR; R.Slope = slope; R.Group = G;
R.PerFly = table(names(:), recH, slope, pf, rz(:, 1), rz(:, 2), 'VariableNames', ...
    {'Fly','RecordHours','Slope','Slope_P_Reversion','RunsZ_Sleep','RunsZ_Wake'});
R.Settings = opts; R.Settings.FlyNames = names; R.Settings.Version = '2026-11-06';
end

% =====================================================================
function V = vt(fl, L, step)
% Variance over start times of the sleep surplus change over each scale.
t0 = fl.Onset(1); tEnd = fl.Onset(end) + fl.Dur(end);
b = [fl.Onset(:); tEnd];                                   % bout boundaries
cs = [0; cumsum(fl.Dur(:) .* fl.IsSleep(:))];              % cumulative sleep at boundaries
pbar = cs(end) / (tEnd - t0);
V = NaN(1, numel(L));
for k = 1:numel(L)
    ts = (t0:step:tEnd - L(k)).';
    if numel(ts) < 10, continue, end
    d = interp1(b, cs, ts + L(k)) - interp1(b, cs, ts) - pbar * L(k);
    V(k) = var(d);
end
end

function z = runs(fl)
z = NaN(1, 2);
types = [true false];
for k = 1:2
    v = fl.Dur(fl.IsSleep == types(k));
    m = median(v);
    a = v(v ~= m) > m;
    if numel(a) < 10, continue, end
    z(k) = 1 + sum(a(2:end) ~= a(1:end-1));
end
end
