function R = ProcessSProxy(flies, opts)
%PROCESSSPROXY  Does a Borbely-type sleep-pressure proxy predict the next bout?
%
%   R = ProcessSProxy(flies)
%   R = ProcessSProxy(flies, TauW=[0.5 1 2 4 8 16], TauS=[0.5 1 2 4 8 16], ...)
%
% Per fly, a pressure proxy S in [0, 1] is run through the fly's own
% record, exactly bout by bout:
%   during wake   S <- 1 - (1 - S)*exp(-duration/TauW)    (pressure builds)
%   during sleep  S <- S*exp(-duration/TauS)              (pressure drains)
% starting at S = 0.5. At each bout's onset S is recorded; bouts starting in
% the first BurnInHours of a contiguous segment are excluded (S still
% depends on the arbitrary start). Then, per grid pair (TauW, TauS):
%   rho_S = Spearman(S at sleep onset, that sleep bout's duration)   homeostat > 0
%   rho_W = Spearman(S at wake onset, that wake bout's duration)     homeostat < 0
%
% GROUP TEST (max-statistic, corrects for scanning the grid): the group
% mean of rho_S, MAXIMISED over the grid, against the same maximum of each
% surrogate index's group means; one-sided. rho_W the same way, minimised.
% Null: RenewalSurrogate (sleep and wake durations permuted within type,
% timeline rebuilt) -- with the same Seed, IDENTICAL surrogates to
% SleepAmountACF. The best grid pair is reported as descriptive: a peak on a
% grid is not an estimate with an interval.
%
% A non-adjacent break (a gap) ends a segment: S restarts at 0.5 and the
% burn-in applies again.
%
% INPUT   flies: struct array (or cell) with Dur (s), IsSleep, Onset (s).
% OPTIONS TauW, TauS (hours; default [0.5 1 2 4 8 16] each); BurnInHours 12;
%         NumSurrogates 1000; Seed 1; FlyNames.
% OUTPUT  struct R
%   RhoS, RhoW        nFly x nW x nS (TauW along dim 2, TauS along dim 3)
%   GroupRhoS, GroupRhoW   nW x nS
%   Group   table: MaxRhoS, BestTauW_S, BestTauS_S, MaxRhoS_P, MinRhoW,
%           BestTauW_W, BestTauS_W, MinRhoW_P, NFlies
%   PerFly  table: Fly, NSleepUsed, NWakeUsed, RhoS_atGroupBest,
%           RhoW_atGroupBest
%   GroupNullMaxS, GroupNullMinW   1 x NumSurrogates
%
% MEASURED CALIBRATION (test_ProcessSProxy, 2026-11-06; 48-fly groups of
% SimulateSleepWake flies): independent bouts (60 groups) rho_S 0.017, rho_W
% 0.017 (conservative); homeostat (30) 1.000 and 1.000, best grid pair at the
% truth (2, 1) h in 9/30; shared slow state with no homeostat (30) rho_S 1.000,
% rho_W 0.367, best (0.5, 0.5) h in 10/30. Calibrated but NOT specific: a
% positive result means a slow variable links successive bouts.
%
% Needs no toolbox.
%
% CHANGE LOG
%   2026-11-06  First version (PRESPEC 2026-11-06 entry D, committed fe307b0).

arguments
    flies
    opts.TauW (1,:) double {mustBePositive} = [0.5 1 2 4 8 16]
    opts.TauS (1,:) double {mustBePositive} = [0.5 1 2 4 8 16]
    opts.BurnInHours (1,1) double {mustBeNonnegative} = 12
    opts.NumSurrogates (1,1) double {mustBeInteger, mustBePositive} = 1000
    opts.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.FlyNames string = string.empty
end
if iscell(flies), flies = [flies{:}]; end
nF = numel(flies);
if isempty(opts.FlyNames), names = "fly" + string(1:nF); else, names = opts.FlyNames(:).'; end
[TW, TS] = ndgrid(opts.TauW * 3600, opts.TauS * 3600);
tw = TW(:).'; ts = TS(:).';
nG = numel(tw);
Rn = opts.NumSurrogates;

rhoS = NaN(nF, nG); rhoW = rhoS; nullS = NaN(nF, nG, Rn); nullW = nullS;
nSu = zeros(nF, 1); nWu = nSu;
for f = 1:nF
    fl = flies(f);
    [rhoS(f, :), rhoW(f, :), nSu(f), nWu(f)] = rhos(fl, tw, ts, opts.BurnInHours * 3600);
    rs = RandStream('mrg32k3a', 'Seed', opts.Seed);
    rs.Substream = f;
    for r = 1:Rn
        [nullS(f, :, r), nullW(f, :, r)] = rhos(RenewalSurrogate(fl, rs), tw, ts, opts.BurnInHours * 3600);
    end
end

gS = mean(rhoS, 1, 'omitnan'); gW = mean(rhoW, 1, 'omitnan');
nMaxS = max(squeeze(mean(nullS, 1, 'omitnan')), [], 1);
nMinW = min(squeeze(mean(nullW, 1, 'omitnan')), [], 1);
[mxS, iS] = max(gS); [mnW, iW] = min(gW);
G = table(mxS, tw(iS)/3600, ts(iS)/3600, (1 + sum(nMaxS >= mxS)) / (Rn + 1), ...
    mnW, tw(iW)/3600, ts(iW)/3600, (1 + sum(nMinW <= mnW)) / (Rn + 1), nF, ...
    'VariableNames', {'MaxRhoS','BestTauW_S','BestTauS_S','MaxRhoS_P', ...
    'MinRhoW','BestTauW_W','BestTauS_W','MinRhoW_P','NFlies'});
sz = [nF, numel(opts.TauW), numel(opts.TauS)];
R.RhoS = reshape(rhoS, sz); R.RhoW = reshape(rhoW, sz);
R.GroupRhoS = reshape(gS, sz(2:3)); R.GroupRhoW = reshape(gW, sz(2:3));
R.Group = G;
R.PerFly = table(names(:), nSu, nWu, rhoS(:, iS), rhoW(:, iW), ...
    'VariableNames', {'Fly','NSleepUsed','NWakeUsed','RhoS_atGroupBest','RhoW_atGroupBest'});
R.GroupNullMaxS = nMaxS; R.GroupNullMinW = nMinW;
R.Settings = opts; R.Settings.FlyNames = names; R.Settings.Version = '2026-11-06';
end

% =====================================================================
function [rS, rW, nS, nW] = rhos(fl, tw, ts, burn)
d = fl.Dur(:); isS = logical(fl.IsSleep(:)); on = fl.Onset(:);
n = numel(d); nG = numel(tw);
Sx = NaN(n, nG);
use = false(n, 1);
S = 0.5 * ones(1, nG);
segStart = on(1);
for i = 1:n
    if i > 1 && on(i) ~= on(i-1) + d(i-1)      % a gap: new segment
        S = 0.5 * ones(1, nG);
        segStart = on(i);
    end
    Sx(i, :) = S;
    use(i) = on(i) - segStart >= burn;
    if isS(i)
        S = S .* exp(-d(i) ./ ts);
    else
        S = 1 - (1 - S) .* exp(-d(i) ./ tw);
    end
end
iS = find(use & isS); iW = find(use & ~isS);
nS = numel(iS); nW = numel(iW);
rS = col_spearman(Sx(iS, :), d(iS));
rW = col_spearman(Sx(iW, :), d(iW));
end

function r = col_spearman(X, y)
% Spearman of every column of X with y. Durations y can tie (integer
% seconds) and get average ranks; S is continuous, so its columns are ranked
% by sorting (an exact tie in S would be broken by order -- it would need two
% bouts with identical histories).
n = size(X, 1);
if n < 10, r = NaN(1, size(X, 2)); return, end
ry = tied_rank(y); ry = ry - mean(ry);
[~, ix] = sort(X, 1);
Rx = zeros(size(X));
Rx(ix + (0:size(X, 2) - 1) * n) = repmat((1:n).', 1, size(X, 2));
Rx = Rx - (n + 1) / 2;
den = sqrt(sum(Rx.^2, 1) * sum(ry.^2));
r = (ry.' * Rx) ./ den;
r(den == 0) = NaN;
end

function r = tied_rank(v)
[~, ~, j] = unique(v);
cnt = accumarray(j, 1);
cum = cumsum(cnt);
avg = cum - (cnt - 1) / 2;
r = avg(j);
end
