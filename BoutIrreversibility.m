function R = BoutIrreversibility(bouts, opts)
%BOUTIRREVERSIBILITY  Model-free time-irreversibility tests on bout sequences.
%
%   R = BoutIrreversibility(bouts)
%   R = BoutIrreversibility(bouts, Lags=[1 2 4 8], NumSurrogates=1000, ...)
%
% For each fly, takes ONE bout type (sleep or wake -- never mixed) in time
% order, forms the differences between bouts k apart,
%
%       D_k(t) = x(t+k) - x(t),      every overlapping pair within a block,
%
% and asks whether their distribution is asymmetric. That is the test of a
% TIME ARROW, and it is model-free: for ANY time-reversible process --
% whatever its bout distribution and whatever form its memory takes --
% (x(t), x(t+k)) has the same distribution as (x(t+k), x(t)), so D_k and
% -D_k have the same distribution and D_k is exactly symmetric. The bout
% distribution itself (exponential, mixture, anything) cannot make D_k
% skewed; only the ORDER can. Contrast a residual in bout space (observed
% minus predicted bout), which is skewed even for a reversible process,
% because the bouts are.
%
% Statistics, per fly and lag (all computed on D_k pooled over the fly's
% blocks; a pair never spans two blocks):
%   Bowley      quartile skewness (Q3 + Q1 - 2*Q2)/(Q3 - Q1) of D_k, with
%               Parzen mid-quantiles so integer-second ties do not bias it
%               (same estimator as GeneralizedHyperbolic_MLE's
%               SkewnessBowleyData). In [-1, 1]; 0 under reversibility.
%               Insensitive to the tails.
%   Moment      moment skewness m3/m2^1.5 of D_k. Tail-sensitive. For lag 1
%               and a stationary sequence, E[D^3] = 3(E[x_t^2 x_t+1] -
%               E[x_t x_t+1^2]), the classical time-reversal statistic.
%               DESCRIPTIVE ONLY -- NO P-VALUE: neither null calibrates it
%               (MEASURED below), because it is dominated by a few extreme
%               bouts.
%   UpFraction  fraction of the nonzero D_k that are positive. 0.5 under
%               reversibility. A sawtooth of many small rises and few large
%               drops gives UpFraction > 0.5 together with NEGATIVE skew.
%               The most powerful of the calibrated tests here (MEASURED).
%   Rho         Spearman rank correlation of (x(t), x(t+k)): the MEMORY at
%               lag k. Time-symmetric, so it is NOT a test of direction;
%               reported because skew and memory are separate questions.
%
% NULLS (both computed by default, each from its own random substream):
%   'iaaft'    Iterative amplitude-adjusted Fourier-transform surrogates
%              (Schreiber & Schmitz 1996), made separately for each block:
%              the SAME bout values, rearranged so the autocorrelation is
%              approximately preserved. They are realisations of a
%              time-reversible process (a Gaussian linear process seen
%              through a monotone transform), so they keep the bout
%              distribution AND the memory and remove only the arrow. Use
%              this p-value for the skew statistics.
%   'shuffle'  Random permutation within each block. Removes memory AND
%              direction. Its p-values for Bowley/UpFraction are
%              valid only if there is no memory -- with memory, the
%              statistic's sampling spread is wider than the shuffle
%              spread, so this null can reject too often (MEASURED below).
%              Its p-value for Rho is the memory test. Rho has no IAAFT
%              p-value: IAAFT keeps the autocorrelation by construction.
%   'strata'   PROFILE-MATCHED permutation, run when the Strata option is
%              given (added 2026-11-05 for light:dark data). Every bout
%              carries a stratum label -- e.g. the 1-h ZT bin of its onset
%              -- and a surrogate permutes bouts only among positions with
%              the same label, across ALL the fly's blocks. Each position
%              therefore keeps its time of day, so the fly's day/night
%              profile (and any trend within a phase that follows ZT) is in
%              the null, and only order BEYOND the profile is tested. Like
%              the shuffle it removes memory; its Rho p-value is the test of
%              memory beyond the profile. Strata with a single bout cannot
%              move (reported as StrataMovable, the share of bouts that can).
%              Why it is needed: in LD the profile alone is irreversible --
%              e.g. night sleep bouts lengthening towards dawn, then a
%              drop at lights-on -- and IAAFT, which does not know the
%              clock, attributes that to the order (MEASURED below).
%              GROUP LEVEL: no group-mean p-value is reported for this
%              null. Permuting within ZT bins cuts the real link between a
%              bout's length and the next bout's onset time (a long night
%              bout pushes the next onset later, where bouts are longer),
%              a per-fly bias too small to see one fly at a time but large
%              once averaged: the group mean over-rejected 0.15-0.18 on
%              10-fly LD-like groups (fixing each bout's ZT in advance
%              removed most of it). The group verdict is instead the COUNT
%              of flies with per-fly UpFraction p < 0.05 against
%              Binomial(NFlies, 0.05): Group.UpFraction_NSig_Strata and
%              UpFraction_CountP_Strata.
%              MEASURED (test_BoutIrreversibility 16-18; LD-like flies: 5
%              days, night sleep bouts rising 20 -> 80 min, latent AR phi
%              0.5, 1-h ZT strata; 200 flies, 199 surrogates):
%                per fly, reversible, lag 1-8: UpFraction 0.035-0.050,
%                  Bowley 0.035-0.085, versus IAAFT 0.070-0.100 and
%                  shuffle 0.090-0.135 for UpFraction (the profile fools
%                  them);
%                per-fly power, irreversible: UpFraction 0.535 at lag 1;
%                count test, 40-fly reversible groups: UpFraction 0.030
%                  (2.06 flies significant on average, 2.00 expected);
%                  the Bowley count would have rejected 0.130, so Bowley
%                  has per-fly p-values only.
%   All nulls are SAMPLE-SIZE MATCHED automatically -- each surrogate has
%   the fly's own blocks, lengths and bout values.
%
% P-VALUES are two-sided Monte Carlo values with the +1 correction,
%   p = min(1, 2*min(1 + #{null >= obs}, 1 + #{null <= obs}) / (Rv + 1)),
% Rv = number of finite null draws. No centring is assumed.
%
% GROUP TEST: the group statistic at each lag is the MEAN over flies of the
% per-fly statistic (each fly weighted equally, not by its pair count). Its
% null takes, for each surrogate index r, the mean over the same flies of
% their r-th surrogate values -- an exact Monte Carlo null for the group
% mean, because surrogates are independent across flies.
%
% INPUT
%   bouts   One of:
%           - a numeric VECTOR (either orientation): one fly, one block;
%           - a numeric MATRIX with >= 2 rows and >= 2 columns: one fly per
%             ROW, bouts in time order along the row, NaN padding allowed
%             only at the END of a row (an interior NaN is an error -- pass
%             blocks instead);
%           - a CELL array with one element per FLY, each element either a
%             numeric vector (one block) or a cell array of numeric vectors
%             (that fly's blocks, each a contiguous stretch of recording).
%           A block is the narrowest window a difference may be formed in:
%           split at recording gaps and at any boundary you do not want a
%           pair to span (ZT window, day, startZT segment).
%
% OPTIONS (name-value)
%   Lags           Positive integers. Default [1 2 4 8].
%   NumSurrogates  Surrogates per fly per null. Default 1000. The smallest
%                  attainable p is 2/(NumSurrogates+1).
%   Null           "both" (default), "iaaft", "shuffle" or "strata". The
%                  skipped null's columns are NaN. "both" also runs the
%                  strata null when Strata is given; "strata" runs only it.
%                  A null that is not run costs nothing; IAAFT is slowest.
%   Strata         Same shape as bouts (vector, matrix, or cell per fly of
%                  vectors/blocks), one label per bout, e.g. floor(ZT) of
%                  its onset. Default [] (no strata null). Its p-values are
%                  columns *_P_Strata, and PerFly gains StrataMovable. Uses
%                  substream 1e6 + f, so adding it never changes the other
%                  nulls' results.
%   Seed           Seeds a private RandStream ('mrg32k3a'); the global
%                  stream is never used. Fly f, null j uses substream
%                  2*(f-1)+j, so results for one fly do not depend on which
%                  null was also run. Default 1.
%   MinPairs       A fly/lag with fewer pairs gets NaN and is left out of
%                  the group test at that lag. Default 20.
%   IAAFTMaxIter   Default 200. Stops earlier when the rank order repeats.
%   FlyNames       String array, one per fly. Default "fly1", "fly2", ...
%   ReturnNull     true to return per-fly null draws (R.PerFlyNull).
%                  Default false.
%   Plot           true for a three-panel figure (light theme pinned).
%
% OUTPUT  struct R with fields
%   PerFly      Table, one row per fly x lag: Fly, Lag, NPairs, NBlocks,
%               Bowley, Bowley_P_IAAFT, Bowley_P_Shuffle, Moment,
%               UpFraction, UpFraction_P_IAAFT, UpFraction_P_Shuffle, Rho,
%               Rho_P_Shuffle.
%   Group       Table, one row per lag: Lag, NFlies, then for each
%               statistic its group mean, NPositive (flies above the null
%               centre -- 0 for skew, 0.5 for UpFraction, 0 for Rho), and
%               P_IAAFT and P_Shuffle (Rho: P_Shuffle only; Moment: none).
%   GroupNull   struct: GroupNull.(null).(stat) is nLag x Rv draws of the
%               group mean; GroupNull.(null).([stat 'Band']) is nLag x 2,
%               the 2.5 and 97.5 percentiles.
%   PerFlyNull  (ReturnNull only) PerFlyNull.(null).(stat): nFly x nLag x R.
%   Settings    the options used, plus Version.
%   Figure      handle when Plot is true.
%
% CAVEATS
%   - D_k values overlap (each bout enters two lag-1 differences), so no
%     iid formula applies to them. That is why every p-value here comes from
%     surrogates of the whole sequence, never from an SE.
%   - A non-significant result means THIS asymmetry is absent at THIS lag;
%     it does not prove reversibility.
%   - IAAFT treats each block as periodic; a block whose two ends differ a
%     lot gets slightly distorted low frequencies. Its null is a Gaussian
%     linear process through a monotone transform, so a rejection says the
%     data are not THAT -- for the skew statistics, whose expectation is 0
%     under ANY reversible process, the IAAFT null only calibrates the
%     spread, which is what the test needs.
%   - Short blocks (a few bouts) give IAAFT surrogates close to shuffles.
%   - Nothing here handles censoring: a bout cut by the block edge enters
%     at its truncated length, as in the rest of Icode.
%
% MEASURED CALIBRATION (test_BoutIrreversibility tests 9-12; 200 simulated
% flies per row, n = 150 bouts, 199 surrogates; rejection rate at nominal
% 0.05, range over lags 1, 2, 4, 8; 200 reps give about +-0.03):
%                                  Bowley         UpFraction     Moment
%   iid exponential      IAAFT    0.035-0.065    0.025-0.055    0.060-0.115
%                        shuffle  0.030-0.050    0.020-0.055    0.020-0.040
%   reversible, memory   IAAFT    0.030-0.055    0.025-0.055    0.075-0.110
%   (latent AR phi 0.5,  shuffle  0.040-0.090    0.040-0.055    0.070-0.220
%    Spearman rho1 ~0.48)
%   POWER, irreversible  IAAFT    lag 1 0.285    lag 1 0.875    (0.935, but
%   (AR phi 0.5, exp.             lag 2 0.080    lag 2 0.265     not valid)
%    innovations)
% So: test with Bowley and UpFraction against IAAFT; UpFraction has by far
% the most power at lag 1 against this alternative; power fades by lag 4
% because the simulated memory does. The shuffle null over-rejects the skew
% statistics once memory is present (Bowley 0.090, Moment 0.220 at lag 1).
% Group test (test 12; 100 groups of 10 flies x 80 bouts, reversible with
% memory, IAAFT): rejection 0.020 (Bowley), 0.040 (UpFraction).
% These are properties of THESE simulated processes: real flies may have
% longer memory or regime switches, which no test here covers.
%
% REFERENCES
%   Weiss G (1975) Time-reversibility of linear stochastic processes.
%     J Appl Prob 12:831-836.
%   Theiler J, Eubank S, Longtin A, Galdrikian B, Farmer JD (1992) Testing
%     for nonlinearity in time series: the method of surrogate data.
%     Physica D 58:77-94.
%   Schreiber T, Schmitz A (1996) Improved surrogate data for nonlinearity
%     tests. Phys Rev Lett 77:635-638.
%   Schreiber T, Schmitz A (2000) Surrogate time series. Physica D
%     142:346-382.
%   Parzen E (2004) Quantiles, conditional quantiles, and mid-quantiles.
%     Statist Sci 19:652-662.
%
% Needs no toolbox.
%
% CHANGE LOG
%   2026-11-05  First version.
%   2026-11-05  Strata option: the profile-matched null for LD data.

arguments
    bouts
    opts.Lags (1,:) double {mustBeInteger, mustBePositive} = [1 2 4 8]
    opts.NumSurrogates (1,1) double {mustBeInteger, mustBePositive} = 1000
    opts.Null (1,1) string {mustBeMember(opts.Null, ["both","iaaft","shuffle","strata"])} = "both"
    opts.Strata = []
    opts.Seed (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.MinPairs (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.IAAFTMaxIter (1,1) double {mustBeInteger, mustBePositive} = 200
    opts.FlyNames string = string.empty
    opts.ReturnNull (1,1) logical = false
    opts.Plot (1,1) logical = false
end

flies = normalise_input(bouts);
nF = numel(flies);
lags = unique(opts.Lags);
nL = numel(lags);
Rn = opts.NumSurrogates;

if isempty(opts.FlyNames)
    names = "fly" + string(1:nF);
else
    names = opts.FlyNames(:).';
    if numel(names) ~= nF
        error('BoutIrreversibility:FlyNames', ...
            'FlyNames has %d entries for %d flies.', numel(names), nF);
    end
end

stats  = ["Bowley","Moment","UpFraction","Rho"];
nulls  = ["iaaft","shuffle","strata"];
hasStrata = ~isempty(opts.Strata);
if opts.Null == "strata" && ~hasStrata
    error('BoutIrreversibility:NoStrata', 'Null="strata" needs the Strata option.');
end
if hasStrata
    strata = normalise_input(opts.Strata);
    if numel(strata) ~= nF
        error('BoutIrreversibility:Strata', 'Strata has %d flies, bouts has %d.', ...
            numel(strata), nF);
    end
    for f = 1:nF
        if numel(strata{f}) ~= numel(flies{f}) || ...
                any(cellfun(@numel, strata{f}) ~= cellfun(@numel, flies{f}))
            error('BoutIrreversibility:Strata', ...
                'Strata of fly %d do not match its blocks bout for bout.', f);
        end
    end
end
useNull = [ismember(opts.Null, ["both","iaaft"]), ismember(opts.Null, ["both","shuffle"]), ...
           hasStrata && ismember(opts.Null, ["both","strata"])];
centre = [0 0 0.5 0];

obs    = NaN(nF, nL, 4);
nPairs = zeros(nF, nL);
nBlk   = zeros(nF, 1);
movable = NaN(nF, 1);
nullDraws = NaN(nF, nL, 4, Rn, 3);

for f = 1:nF
    blocks = flies{f};
    nBlk(f) = numel(blocks);
    [S, np] = seq_stats(blocks, lags);
    nPairs(f, :) = np;
    ok = np >= opts.MinPairs;
    S(~ok, :) = NaN;
    obs(f, :, :) = reshape(S, 1, nL, 4);
    if hasStrata
        % Positions sharing a stratum, over all the fly's blocks: the
        % profile-matched surrogate permutes bouts only within these sets,
        % so every position keeps its stratum (e.g. its time of day).
        xAll = vertcat(blocks{:});
        [~, ~, g] = unique(vertcat(strata{f}{:}));
        groups = accumarray(g, (1:numel(g)).', [], @(v) {v});
        groups = groups(cellfun(@numel, groups) > 1);
        movable(f) = sum(cellfun(@numel, groups)) / numel(xAll);
        blen = cellfun(@numel, blocks);
    end
    for j = 1:3
        if ~useNull(j), continue, end
        rs = RandStream('mrg32k3a', 'Seed', opts.Seed);
        if j < 3
            rs.Substream = 2*(f-1) + j;
        else
            rs.Substream = 1e6 + f;    % kept apart from the first two nulls'
        end
        for r = 1:Rn
            if j == 3
                y = xAll;
                for k = 1:numel(groups)
                    ix = groups{k};
                    y(ix) = xAll(ix(randperm(rs, numel(ix))));
                end
                sb = mat2cell(y, blen, 1).';
            else
                sb = cell(size(blocks));
                for b = 1:numel(blocks)
                    if j == 1
                        sb{b} = iaaft(blocks{b}, rs, opts.IAAFTMaxIter);
                    else
                        xb = blocks{b};
                        sb{b} = xb(randperm(rs, numel(xb)));
                    end
                end
            end
            Sr = seq_stats(sb, lags);
            Sr(~ok, :) = NaN;
            nullDraws(f, :, :, r, j) = reshape(Sr, 1, nL, 4);
        end
    end
end
% ---- per-fly p-values ----
% IAAFT keeps the autocorrelation, so it is no null for Rho: no p-value,
% but its Rho draws are kept (ReturnNull) as a check that it did.
P = NaN(nF, nL, 4, 3);
for j = 1:3
    if ~useNull(j), continue, end
    for s = 1:4
        if ~hasP(s, j), continue, end
        for f = 1:nF
            for l = 1:nL
                P(f, l, s, j) = mc_p(obs(f, l, s), squeeze(nullDraws(f, l, s, :, j)));
            end
        end
    end
end

% ---- per-fly table ----
[FL, LL] = ndgrid(1:nF, 1:nL);
FL = FL.'; LL = LL.';
FL = FL(:); LL = LL(:);
idx = sub2ind([nF nL], FL, LL);
col = @(v) reshape(v, [], 1);
T = table(col(names(FL)), col(lags(LL)), col(nPairs(idx)), col(nBlk(FL)), ...
    'VariableNames', {'Fly','Lag','NPairs','NBlocks'});
for s = 1:4
    o = obs(:, :, s);
    T.(stats(s)) = col(o(idx));
    for j = 1:3
        if ~hasP(s, j) || (j == 3 && ~hasStrata), continue, end
        p = P(:, :, s, j);
        T.(stats(s) + "_P_" + label(nulls(j))) = col(p(idx));
    end
end
if hasStrata
    T.StrataMovable = col(movable(FL));
end
R.PerFly = T;

% ---- group test ----
G = table(lags.', 'VariableNames', {'Lag'});
G.NFlies = sum(isfinite(obs(:, :, 1)), 1).';
GroupNull = struct();
for j = 1:3
    if j == 3 && ~hasStrata, continue, end
    for s = 1:4
        draws = NaN(nL, Rn);
        gp = NaN(nL, 1);
        gm = NaN(nL, 1);
        npos = NaN(nL, 1);
        for l = 1:nL
            o = obs(:, l, s);
            in = isfinite(o);
            if ~any(in), continue, end
            gm(l) = mean(o(in));
            npos(l) = sum(o(in) > centre(s));
            if useNull(j) && hasP(s, j)
                nd = squeeze(nullDraws(in, l, s, :, j));
                if sum(in) == 1, nd = nd(:).'; end
                draws(l, :) = mean(nd, 1, 'omitnan');
                gp(l) = mc_p(gm(l), draws(l, :).');
            end
        end
        if j == 1
            G.(stats(s)) = gm;
            G.(stats(s) + "_NPositive") = npos;
        end
        if hasP(s, j) && j ~= 3   % no group p for strata: biased (see header)
            G.(stats(s) + "_P_" + label(nulls(j))) = gp;
        end
        if s == 3 && j == 3 && useNull(3)
            % Group verdict for the strata null: COUNT of flies with
            % per-fly p < 0.05, against Binomial(NFlies, 0.05), upper tail.
            % UpFraction only: measured calibrated (0.030 on 40-fly LD-like
            % groups); the Bowley count over-rejected (0.130, per-fly rate
            % about 0.068), so Bowley gets no count verdict.
            pf = P(:, :, s, 3);
            nsig = sum(pf < 0.05, 1).';
            ntest = sum(isfinite(pf), 1).';
            G.(stats(s) + "_NSig_Strata") = nsig;
            G.(stats(s) + "_CountP_Strata") = arrayfun(@(k, n) binom_upper(k, n, 0.05), nsig, ntest);
        end
        GroupNull.(nulls(j)).(stats(s)) = draws;
        GroupNull.(nulls(j)).(stats(s) + "Band") = band(draws);
    end
end
R.Group = G;
R.GroupNull = GroupNull;
if opts.ReturnNull
    for j = 1:2 + hasStrata
        for s = 1:4
            R.PerFlyNull.(nulls(j)).(stats(s)) = squeeze(nullDraws(:, :, s, :, j));
        end
    end
end
R.Settings = opts;
R.Settings.Lags = lags;
R.Settings.FlyNames = names;
R.Settings.Version = '2026-11-05';
if opts.Plot
    R.Figure = plot_result(R, obs, lags, centre);
end
end

% =====================================================================
function flies = normalise_input(bouts)
if isnumeric(bouts)
    if isvector(bouts)
        flies = {{double(bouts(:))}};
    else
        flies = cell(size(bouts, 1), 1);
        for i = 1:size(bouts, 1)
            flies{i} = {strip_pad(double(bouts(i, :)), i)};
        end
    end
elseif iscell(bouts)
    flies = cell(numel(bouts), 1);
    for i = 1:numel(bouts)
        e = bouts{i};
        if isnumeric(e)
            flies{i} = {check_block(double(e(:)), i)};
        elseif iscell(e)
            flies{i} = cellfun(@(b) check_block(double(b(:)), i), e(:).', ...
                'UniformOutput', false);
        else
            error('BoutIrreversibility:Input', ...
                'Fly %d is neither numeric nor a cell array of blocks.', i);
        end
    end
else
    error('BoutIrreversibility:Input', ...
        'bouts must be numeric or a cell array (see help).');
end
end

function x = strip_pad(row, i)
last = find(~isnan(row), 1, 'last');
if isempty(last), last = 0; end
x = check_block(row(1:last).', i);
end

function x = check_block(x, i)
if any(isnan(x))
    error('BoutIrreversibility:InteriorNaN', ...
        ['Fly %d has a NaN inside a block. Split the sequence into ' ...
         'blocks at the gap instead.'], i);
end
if any(~isfinite(x))
    error('BoutIrreversibility:NonFinite', 'Fly %d has a non-finite bout.', i);
end
end

% =====================================================================
function [S, np] = seq_stats(blocks, lags)
% S: nLag x 4 [Bowley Moment UpFraction Rho]; np: pairs per lag.
nL = numel(lags);
S = NaN(nL, 4);
np = zeros(1, nL);
for l = 1:nL
    k = lags(l);
    a = cell(numel(blocks), 1);
    c = a;
    for b = 1:numel(blocks)
        x = blocks{b};
        if numel(x) > k
            a{b} = x(1:end-k);
            c{b} = x(1+k:end);
        end
    end
    a = vertcat(a{:});
    c = vertcat(c{:});
    np(l) = numel(a);
    if np(l) < 2, continue, end
    d = c - a;
    S(l, 1) = bowley_mid(d);
    dc = d - mean(d);
    m2 = mean(dc.^2);
    if m2 > 0, S(l, 2) = mean(dc.^3) / m2^1.5; end
    nz = d ~= 0;
    if any(nz), S(l, 3) = mean(d(nz) > 0); end
    S(l, 4) = spearman(a, c);
end
end

function b = bowley_mid(x)
% Bowley's quartile skewness with Parzen mid-quantiles; identical to the
% estimator in GeneralizedHyperbolic_MLE (SkewnessBowleyData).
x = sort(x(:));
n = numel(x);
[u, ~, j] = unique(x);
b = NaN;
if numel(u) < 2, return, end
c  = accumarray(j, 1);
Fm = (cumsum(c) - c/2) / n;
P  = [0.25 0.5 0.75];
q  = interp1(Fm, u, P, 'linear');
q(P < Fm(1))   = u(1);
q(P > Fm(end)) = u(end);
d = q(3) - q(1);
if d > 0, b = (q(3) + q(1) - 2*q(2)) / d; end
end

function r = spearman(a, c)
ra = tied_rank(a);
rc = tied_rank(c);
ra = ra - mean(ra);
rc = rc - mean(rc);
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

% =====================================================================
function y = iaaft(x, rs, maxIter)
% Iterative amplitude-adjusted Fourier-transform surrogate
% (Schreiber & Schmitz 1996): same values as x, approximately the same
% power spectrum. Ends on the rank-ordering step, so the values are exact.
x = x(:);
n = numel(x);
if n < 3
    y = x(randperm(rs, n));
    return
end
xs = sort(x);
A = abs(fft(x));
y = x(randperm(rs, n));
[~, prev] = sort(y);
for it = 1:maxIter
    Y = fft(y);
    y = real(ifft(A .* exp(1i * angle(Y))));
    [~, ix] = sort(y);
    y(ix) = xs;
    if isequal(ix, prev), break, end
    prev = ix;
end
end

% =====================================================================
function p = mc_p(o, null)
null = null(isfinite(null));
if ~isfinite(o) || isempty(null)
    p = NaN;
    return
end
Rv = numel(null);
p = min(1, 2 * min(1 + sum(null >= o), 1 + sum(null <= o)) / (Rv + 1));
end

function b = band(draws)
b = NaN(size(draws, 1), 2);
for l = 1:size(draws, 1)
    d = sort(draws(l, isfinite(draws(l, :))));
    if isempty(d), continue, end
    if isscalar(d), b(l, :) = d; continue, end
    b(l, :) = interp1(linspace(0, 1, numel(d)), d, [0.025 0.975], 'linear');
end
end

function s = label(nullName)
s = containers.Map({'iaaft','shuffle','strata'}, {"IAAFT","Shuffle","Strata"});
s = s(char(nullName));
end

% =====================================================================
function fh = plot_result(R, obs, lags, centre)
fh = figure('Name', 'BoutIrreversibility', 'Color', 'w');
try theme(fh, 'light'); catch, end %#ok<CTCH>
titles = ["Skewness of the change between bouts k apart (Bowley)", ...
          "Fraction of changes that are increases", ...
          "Memory: rank correlation of bouts k apart"];
which = [1 3 4];
names = ["Bowley","UpFraction","Rho"];
tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact');
for i = 1:3
    ax = nexttile(tl);
    hold(ax, 'on');
    s = which(i);
    nullsToDraw = ["iaaft","shuffle"];
    cols = [0.35 0.55 0.85; 0.85 0.55 0.35];
    for j = 1:2
        if s == 4 && j == 1, continue, end
        b = R.GroupNull.(nullsToDraw(j)).(names(i) + "Band");
        if all(isnan(b(:))), continue, end
        fill(ax, [lags fliplr(lags)], [b(:, 1).' fliplr(b(:, 2).')], cols(j, :), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'DisplayName', "95% of group means, " + nullsToDraw(j) + " null");
    end
    for f = 1:size(obs, 1)
        plot(ax, lags, obs(f, :, s), '-', 'Color', [0.75 0.75 0.75], ...
            'HandleVisibility', 'off');
    end
    plot(ax, lags, R.Group.(names(i)), 'k-o', 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'k', 'DisplayName', 'mean over flies');
    yline(ax, centre(s), ':', 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log', 'XTick', lags);
    xlabel(ax, 'Lag k (bouts)');
    title(ax, titles(i), 'FontWeight', 'normal');
    if i == 1, legend(ax, 'Location', 'best'); end
    box(ax, 'on');
end
end

function tf = hasP(s, j)
% Which statistic/null pairs get a p-value. s: 1 Bowley, 2 Moment,
% 3 UpFraction, 4 Rho; j: 1 IAAFT, 2 shuffle.
%   Moment: NONE. Measured miscalibrated with both nulls (test suite,
%     200 simulated flies, n = 150): IAAFT rejected 0.095-0.115 at nominal
%     0.05 even on iid bouts; shuffle rejected 0.22 at lag 1 with memory.
%   Rho: shuffle only (IAAFT keeps the autocorrelation by construction).
tf = s ~= 2 && ~(s == 4 && j == 1);
end

function p = binom_upper(k, n, q)
% P(X >= k) for X ~ Binomial(n, q), toolbox-free.
if n == 0, p = NaN; return, end
j = k:n;
lp = gammaln(n + 1) - gammaln(j + 1) - gammaln(n - j + 1) + j*log(q) + (n - j)*log1p(-q);
p = min(1, sum(exp(lp)));
end
