function r = loopPersistenceFly(y, ph, options)
%LOOPPERSISTENCEFLY  Persistence of loop direction in one fly against clock-matched and 12-h chain nulls (PRESPEC 26; sqw26.py fly()).
%
%   r = loopPersistenceFly(y, ph)
%   r = loopPersistenceFly(y, ph, NSim=500, Window=240, RandomSeed=2600, NShuffle=2000)
%
%   y    (T x 1) minute symbols in the port's coding (1 = Q, 2 = W, 3 = S, 0 = missing), whole days from ZT0
%        (or from the record start in DD)
%   ph   (T x 1) phase of each minute (MINUTEPHASES)
%
%   Loops: completed cycles in order (loop erasure as COMPLETEDLOOPS, without timing), each with a
%   segment number (a missing minute starts a new segment). S = fraction of consecutive loop pairs IN
%   THE SAME SEGMENT with the same direction; runs = lengths of same-direction runs.
%   Nulls, NSim records each, same length, starting in the first scored state, missing minutes blanked:
%     clock  one transition matrix per clock window of Window minutes (6 per day at 240), pooled over
%            days, rows normalised; an empty row takes the 12-h phase matrix of the window's first minute
%     phase  one matrix per 12-h phase (DD: one); an empty row is uniform
%   E = S - mean null S (finite records). Shuffle: mean S over NShuffle permutations of the observed
%   directions with the segments kept (descriptive).
%   r    struct: nLoops, nPairs, fracFwd, S_obs, runs_obs, loopT (completion minutes), loopDir,
%        S_null_clock, S_null_sd_clock, E_clock, runs_null_clock (first 100 records), and the same
%        for phase; S_shuffle. Monte Carlo: E and S_null equal the Python's only within simulation error;
%        nLoops, nPairs, S_obs and runs_obs are exact.
%
%   Example (from the Icode\sqw folder):
%       [y, ph] = sqwFlyInput("D:\flies\per0DD\fly1.txt", false);
%       r = loopPersistenceFly(y, ph); [r.S_obs r.E_clock]
%
%   See also LOOPPERSISTENCETESTS, EXPORTLOOPPERSISTENCEFIGURE, COMPLETEDLOOPS.

arguments
    y (:,1) double
    ph (:,1) double
    options.NSim (1,1) double {mustBePositive, mustBeInteger} = 500
    options.Window (1,1) double {mustBePositive, mustBeInteger} = 240
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 2600
    options.NShuffle (1,1) double {mustBeNonnegative, mustBeInteger} = 2000
end
rs = RandStream('twister', 'Seed', options.RandomSeed);
y(y == 4) = 0; T = numel(y); miss = y == 0;
[seg, dir, tl] = loopDirs(y); [S, np] = Sstat(seg, dir);
r = struct('nLoops', numel(dir), 'nPairs', np, 'fracFwd', mean(dir == 1), 'S_obs', S, 'runs_obs', runLengths(seg, dir), ...
    'loopT', tl, 'loopDir', dir);
nP = max(ph); Pph = chainMatrices(y, ph, nP, []);
wpos = floor(mod((0:T - 1)', 1440) / options.Window) + 1; nw = 1440 / options.Window;
fb = zeros(3, 3, nw); for m = 1:nw, fb(:, :, m) = Pph(:, :, ph(find(wpos == m, 1))); end
Pclk = chainMatrices(y, wpos, nw, fb);
y0 = y(find(~miss, 1)); if isempty(y0), y0 = 1; end
nulls = {"clock", Pclk, wpos; "phase", Pph, ph};
for k = 1:2
    X = simulateChain(nulls{k, 2}, nulls{k, 3}, y0, T, options.NSim, rs); X(:, miss) = 0;
    Ss = NaN(options.NSim, 1); rl = [];
    for s = 1:options.NSim
        [sg, dr] = loopDirs(X(s, :)'); Ss(s) = Sstat(sg, dr);
        if s <= 100, rl = [rl; runLengths(sg, dr)]; end %#ok<AGROW>
    end
    nm = nulls{k, 1}; mu = mean(Ss, 'omitnan');
    r.("S_null_" + nm) = mu; r.("S_null_sd_" + nm) = std(Ss(isfinite(Ss)), 1); r.("E_" + nm) = S - mu; r.("runs_null_" + nm) = rl;
end
r.S_shuffle = NaN;
if np > 0
    sh = NaN(options.NShuffle, 1);
    for s = 1:options.NShuffle, sh(s) = Sstat(seg, dir(randperm(rs, numel(dir)))); end
    r.S_shuffle = mean(sh, 'omitnan');
end
end

function [seg, dir, tl] = loopDirs(y)
% loops in completion order on the run-compressed symbols; a missing run starts a new segment
ch = [1; find(diff(y) ~= 0) + 1]; q = y(ch);
seg = zeros(0, 1); dir = seg; tl = seg; st = zeros(3, 1); n = 0; sg = 0;
for k = 1:numel(q)
    x = q(k);
    if x == 0, n = 0; sg = sg + 1; continue, end
    j = find(st(1:n) == x, 1);
    if isempty(j), n = n + 1; st(n) = x; continue, end
    if n - j + 1 == 3
        lp = st(j:n)'; f = isequal(lp, [2 1 3]) || isequal(lp, [1 3 2]) || isequal(lp, [3 2 1]);
        seg(end + 1, 1) = sg; dir(end + 1, 1) = 2 * f - 1; tl(end + 1, 1) = ch(k); %#ok<AGROW>
    end
    n = j;
end
end

function [S, n] = Sstat(seg, dir)
same = seg(1:end - 1) == seg(2:end); n = nnz(same);
if n == 0, S = NaN; else, S = nnz(same & dir(1:end - 1) == dir(2:end)) / n; end
end

function r = runLengths(seg, dir)
if isempty(dir), r = zeros(0, 1); return, end
b = [true; seg(2:end) ~= seg(1:end - 1) | dir(2:end) ~= dir(1:end - 1)];
r = diff([find(b); numel(dir) + 1]);
end

function P = chainMatrices(y, idx, nm, fb)
% transition matrix per index of the ARRIVAL minute, between consecutive scored minutes
C = zeros(3, 3, nm); a = y(1:end - 1); b = y(2:end); m = idx(2:end); ok = a > 0 & b > 0;
C = C + reshape(accumarray(sub2ind([3 3 nm], a(ok), b(ok), m(ok)), 1, [9 * nm 1]), 3, 3, nm);
P = zeros(3, 3, nm);
for k = 1:nm
    for s = 1:3
        rr = sum(C(s, :, k));
        if rr > 0, P(s, :, k) = C(s, :, k) / rr; elseif isempty(fb), P(s, :, k) = 1 / 3; else, P(s, :, k) = fb(s, :, k); end
    end
end
end

function X = simulateChain(P, idx, y0, T, R, rs)
Cc = cumsum(P, 2); Cc(:, 3, :) = 1; X = zeros(R, T); cur = repmat(y0, R, 1); X(:, 1) = cur;
for t = 2:T
    C = Cc(:, :, idx(t)); u = rand(rs, R, 1);
    cur = 1 + (u > C(cur, 1)) + (u > C(cur, 2)); X(:, t) = cur;
end
end
