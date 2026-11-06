function R = BoutCrossType(flies, opts)
%BOUTCROSSTYPE  Sleep against the wake bout before it and after it (homeostat test).
%
%   R = BoutCrossType(flies)
%   R = BoutCrossType(flies, Null="strata", NumSurrogates=1000, ...)
%
% A homeostat predicts that a long wake bout builds pressure and so is
% followed by a long sleep bout. Per fly, over ADJACENT triples (wake
% w_pre, sleep s, wake w_post) -- adjacent meaning each bout starts exactly
% where the previous one ends -- this computes
%   r_pre  = Spearman(s, w_pre)       homeostat: > 0
%   r_post = Spearman(s, w_post)
%   Delta  = r_pre - r_post           reversing time swaps r_pre and r_post,
%                                     so Delta is a time-arrow statistic
%   r_pre2 = Spearman(s, wake two bouts back)
% and the mirror for wake bouts: q_pre = Spearman(w, s_pre), q_post =
% Spearman(w, s_post), DeltaW = q_pre - q_post.
%
% NULLS (pairing broken, each sequence's own structure kept)
%   "shift"   (default; for data with no time-of-day profile) the fly's
%             wake sequence is circularly shifted against its sleep
%             sequence (sleep against wake for the mirror). Per-fly
%             p-values use EVERY shift k with 2 <= k <= m-2 (m = number of
%             wake bouts). The group null uses NumSurrogates random shifts
%             per fly (k drawn uniformly from that range), so the r-th draws
%             of all flies form one group-mean null draw.
%   "strata"  (data with a profile) the fly's wake bouts are permuted among
%             wake positions with the same Strata label (e.g. 1-h ZT bin of
%             onset), across days, sleep fixed; NumSurrogates draws. The
%             group verdict is then the COUNT of flies with per-fly p < 0.05
%             against Binomial(n, 0.05): no group-mean p is given, because
%             permuting a wake bout cuts its link with the next sleep's
%             onset time (the timeline-coupling bias measured for the strata
%             null in BoutIrreversibility).
% P-values: one-sided (greater) for r_pre, r_pre2, Delta, q_pre, DeltaW;
% two-sided for r_post and q_post. Monte Carlo with the +1 correction.
%
% INPUT
%   flies  struct array (or cell array of structs), one per fly, with
%          column vectors Dur (s), IsSleep (logical), Onset (s); Strata
%          (labels) for Null="strata". Bouts in time order.
% OPTIONS
%   Null           "shift" (default) or "strata".
%   NumSurrogates  Default 1000.
%   Seed           Private RandStream ('mrg32k3a', substream = fly index).
%   MinTriples     A fly with fewer sleep triples gets NaN. Default 20.
%   FlyNames       Default "fly1", ...
% OUTPUT  struct R
%   PerFly  table: Fly, NTriples, NWakeTriples, NNonAdjacent, r_pre,
%           r_pre_P, r_post, r_post_P, Delta, Delta_P, r_pre2, r_pre2_P,
%           q_pre, q_pre_P, q_post, q_post_P, DeltaW, DeltaW_P.
%   Group   table, one row per statistic: Stat, NFlies, Mean, MeanP (shift
%           only; one- or two-sided as above), NSig (per-fly p < 0.05),
%           CountP (Binomial(NFlies, 0.05) upper tail).
%   Settings
%
% MEASURED CALIBRATION (test_BoutCrossType, 2026-11-06; SimulateSleepWake flies;
% rejection at nominal 0.05). Per fly (200), independent bouts: shift 0.025-0.075,
% strata 0.020-0.075 over the seven statistics. 48-fly groups (100), shift, group
% mean: r_pre 0.050, Delta 0.040 (r_pre2 0.100). Shared slow state, no homeostat,
% per fly: r_pre 0.375 but Delta 0.020. LD-like profile, strata, per fly: r_pre
% 0.080, Delta 0.090; 40-fly COUNT: r_pre 0.060 but Delta 0.120, DeltaW 0.150 (use
% the Delta count descriptively only). Homeostat (Tw 2 h, Ts 1 h, gain 2), 48-fly
% groups (50): r_pre 1.000, Delta 0.980. So with the shift null, Delta separates
% a homeostat from a shared slow state; r_pre alone does not.
%
% Needs no toolbox.
%
% CHANGE LOG
%   2026-11-06  First version (PRESPEC 2026-11-06 entry A, committed fe307b0).

arguments
    flies
    opts.Null (1,1) string {mustBeMember(opts.Null, ["shift","strata"])} = "shift"
    opts.NumSurrogates (1,1) double {mustBeInteger, mustBePositive} = 1000
    opts.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.MinTriples (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.FlyNames string = string.empty
end
if iscell(flies), flies = [flies{:}]; end
nF = numel(flies);
if isempty(opts.FlyNames), names = "fly" + string(1:nF); else, names = opts.FlyNames(:).'; end
Rn = opts.NumSurrogates;
stats = ["r_pre","r_post","Delta","r_pre2","q_pre","q_post","DeltaW"];
oneSided = [true false true true true false true];
nS = numel(stats);
obs = NaN(nF, nS); P = NaN(nF, nS); gdraw = NaN(nF, nS, Rn);
nT = zeros(nF, 1); nTW = nT; nNA = nT;

for f = 1:nF
    fl = flies(f);
    d = fl.Dur(:); isS = logical(fl.IsSleep(:)); on = fl.Onset(:);
    n = numel(d);
    adj = [on(2:end) == on(1:end-1) + d(1:end-1); false];   % adj(i): i -> i+1
    nNA(f) = sum(~adj(1:end-1));
    % positions of each type in order, and each bout's index within its type
    posW = find(~isS); posS = find(isS);
    idxIn = zeros(n, 1); idxIn(posW) = 1:numel(posW); idxIn(posS) = 1:numel(posS);
    % sleep triples: i sleep, i-1 and i+1 wake, both links adjacent
    iS = find(isS & [false; adj(1:end-1)] & adj);
    iS = iS(iS > 1 & iS < n);
    iS = iS(~isS(iS - 1) & ~isS(iS + 1));
    % wake two back: i-3 wake with i-3 -> i-2 -> i-1 adjacent
    ok2 = iS > 3;
    ok2(ok2) = adj(iS(ok2) - 3) & adj(iS(ok2) - 2) & ~isS(iS(ok2) - 3);
    iW = find(~isS & [false; adj(1:end-1)] & adj);
    iW = iW(iW > 1 & iW < n);
    iW = iW(isS(iW - 1) & isS(iW + 1));
    nT(f) = numel(iS); nTW(f) = numel(iW);
    if nT(f) < opts.MinTriples, continue, end
    % index (within the wake list / sleep list) of each partner
    jWpre = idxIn(iS - 1); jWpost = idxIn(iS + 1);
    jW2 = NaN(size(iS)); jW2(ok2) = idxIn(iS(ok2) - 3);
    jSpre = idxIn(iW - 1); jSpost = idxIn(iW + 1);
    s = d(iS); w = d(iW);
    Wd = d(posW); Sd = d(posS);
    obs(f, :) = compute(s, w, Wd, Sd, jWpre, jWpost, jW2, jSpre, jSpost);

    rs = RandStream('mrg32k3a', 'Seed', opts.Seed);
    rs.Substream = f;
    mW = numel(Wd); mS = numel(Sd);
    if opts.Null == "shift"
        % per-fly p from every shift
        kk = 2:mW - 2;
        null = NaN(numel(kk), nS);
        for a = 1:numel(kk)
            k = kk(a);
            kS = min(k, mS - 2);
            null(a, :) = compute(s, w, circshift(Wd, k), circshift(Sd, kS), ...
                jWpre, jWpost, jW2, jSpre, jSpost);
        end
        P(f, :) = pvals(obs(f, :), null, oneSided);
        % group draws: random shifts
        for r = 1:Rn
            k = randi(rs, [2, mW - 2]);
            kS = min(k, mS - 2);
            gdraw(f, :, r) = compute(s, w, circshift(Wd, k), circshift(Sd, kS), ...
                jWpre, jWpost, jW2, jSpre, jSpost);
        end
    else
        st = fl.Strata(:);
        gW = groups(st(posW)); gS = groups(st(posS));
        null = NaN(Rn, nS);
        for r = 1:Rn
            Wp = Wd; Sp = Sd;
            for g = 1:numel(gW), ix = gW{g}; Wp(ix) = Wd(ix(randperm(rs, numel(ix)))); end
            for g = 1:numel(gS), ix = gS{g}; Sp(ix) = Sd(ix(randperm(rs, numel(ix)))); end
            null(r, :) = compute(s, w, Wp, Sp, jWpre, jWpost, jW2, jSpre, jSpost);
        end
        P(f, :) = pvals(obs(f, :), null, oneSided);
        gdraw(f, :, :) = permute(null, [3 2 1]);
    end
end

T = table(names(:), nT, nTW, nNA, 'VariableNames', {'Fly','NTriples','NWakeTriples','NNonAdjacent'});
for k = 1:nS
    T.(stats(k)) = obs(:, k);
    T.(stats(k) + "_P") = P(:, k);
end
R.PerFly = T;

G = table(stats(:), 'VariableNames', {'Stat'});
G.NFlies = zeros(nS, 1); G.Mean = NaN(nS, 1); G.MeanP = NaN(nS, 1);
G.NSig = zeros(nS, 1); G.CountP = NaN(nS, 1);
for k = 1:nS
    in = isfinite(obs(:, k));
    G.NFlies(k) = sum(in);
    if ~any(in), continue, end
    G.Mean(k) = mean(obs(in, k));
    if opts.Null == "shift"
        nd = squeeze(mean(gdraw(in, k, :), 1, 'omitnan'));
        G.MeanP(k) = pvals(G.Mean(k), nd(:), oneSided(k));
    end
    G.NSig(k) = sum(P(in, k) < 0.05);
    G.CountP(k) = binom_upper(G.NSig(k), sum(isfinite(P(in, k))), 0.05);
end
R.Group = G;
R.Settings = opts;
R.Settings.FlyNames = names;
R.Settings.Version = '2026-11-06';
end

% =====================================================================
function v = compute(s, w, Wd, Sd, jWpre, jWpost, jW2, jSpre, jSpost)
rpre = spearman(s, Wd(jWpre));
rpost = spearman(s, Wd(jWpost));
ok = isfinite(jW2);
if sum(ok) >= 10, rpre2 = spearman(s(ok), Wd(jW2(ok))); else, rpre2 = NaN; end
if numel(w) >= 10
    qpre = spearman(w, Sd(jSpre)); qpost = spearman(w, Sd(jSpost));
else
    qpre = NaN; qpost = NaN;
end
v = [rpre, rpost, rpre - rpost, rpre2, qpre, qpost, qpre - qpost];
end

function p = pvals(o, null, oneSided)
p = NaN(size(o));
for k = 1:numel(o)
    nk = null(:, k); nk = nk(isfinite(nk));
    if ~isfinite(o(k)) || isempty(nk), continue, end
    if oneSided(k)
        p(k) = (1 + sum(nk >= o(k))) / (numel(nk) + 1);
    else
        p(k) = min(1, 2 * min(1 + sum(nk >= o(k)), 1 + sum(nk <= o(k))) / (numel(nk) + 1));
    end
end
end

function g = groups(labels)
[~, ~, j] = unique(labels);
g = accumarray(j, (1:numel(j)).', [], @(v) {v});
g = g(cellfun(@numel, g) > 1);
end

function r = spearman(a, c)
ra = tied_rank(a); rc = tied_rank(c);
ra = ra - mean(ra); rc = rc - mean(rc);
den = sqrt(sum(ra.^2) * sum(rc.^2));
if den > 0, r = sum(ra .* rc) / den; else, r = NaN; end
end

function r = tied_rank(v)
[~, ~, j] = unique(v);
cnt = accumarray(j, 1);
cum = cumsum(cnt);
avg = cum - (cnt - 1) / 2;
r = avg(j);
end

function p = binom_upper(k, n, q)
if n == 0, p = NaN; return, end
j = k:n;
lp = gammaln(n + 1) - gammaln(j + 1) - gammaln(n - j + 1) + j*log(q) + (n - j)*log1p(-q);
p = min(1, sum(exp(lp)));
end
