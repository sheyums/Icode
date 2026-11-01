function I = influence_20r(W, xmin, label)
% INFLUENCE_20R  The 2026-09-15 per-fly removal design, unchanged ("per0 DD
% wake: does any single fly carry the hazard rise?"), applied to a fly x bout
% matrix W. Run BY NAME from the Icode folder.
%
% For each fly i: RiseRatio with fly i removed, against 199 removals of n_i
% bouts drawn uniformly without replacement from all bouts (rng(1000 + i,
% 'twister'), drawn in sequence). NumBins {12,16,20,24}, Spacing "log",
% MinAtRisk 10, Alpha 0.05, no null replicates, no bootstrap. Fly i is
% influential at a count if its RiseRatio is below the 2.5th percentile
% (prctile) of its 199 removals, and overall at 3 or more counts. If any fly is
% influential: remove all of them, refit the K = 2 null and rerun the rise test
% (rise_test_20r); "carried by a few flies" if RiseRatio <= that run's null
% median at 3 or more counts.

assert(isequal(fileparts(which('BoutHazard')), pwd), ...
    'Run from the Icode folder: BoutHazard resolves to %s', which('BoutHazard'));

nb = [12 16 20 24];
keepB = isfinite(W) & W >= xmin;
nFly = size(W, 1);
fly = repmat((1:nFly)', 1, size(W, 2));
x = W(keepB); f = fly(keepB);
N = numel(x);
rr = @(v, k) BoutHazard(v, xmin, SamplingInterval=1, NumBins=k, Spacing="log", ...
    MinAtRisk=10, Alpha=0.05).RiseRatio;

full = arrayfun(@(k) rr(x, k), nb);
ni = accumarray(f, 1, [nFly 1]);
without = nan(nFly, 4); p025 = nan(nFly, 4); med = nan(nFly, 4);
for i = 1:nFly
    xi = x(f ~= i);
    without(i, :) = arrayfun(@(k) rr(xi, k), nb);
    rng(1000 + i, 'twister');
    rnd = nan(199, 4);
    for b = 1:199
        drop = randperm(N, ni(i));
        xr = x; xr(drop) = [];
        rnd(b, :) = arrayfun(@(k) rr(xr, k), nb);
    end
    p025(i, :) = prctile(rnd, 2.5);
    med(i, :) = median(rnd);
end
below = without < p025;
infl = sum(below, 2) >= 3;

I = struct('label', label, 'n', N, 'nPerFly', ni, 'NumBins', nb, 'full', full, ...
    'without', without, 'p025', p025, 'median', med, 'below', below, 'influential', infl);

fprintf('%s: per-fly influence, n = %d, %d flies\n', label, N, nFly);
fprintf('  full RiseRatio %.3f %.3f %.3f %.3f\n', full);
fprintf('  fly  n_i   without fly i (12/16/20/24)      2.5th pct of random removals    below\n');
for i = 1:nFly
    fprintf('  %3d %4d   %.3f %.3f %.3f %.3f   %.3f %.3f %.3f %.3f   %d%d%d%d\n', ...
        i, ni(i), without(i, :), p025(i, :), below(i, :));
end
fprintf('  influential flies (>= 3 counts): %s\n', mat2str(find(infl)'));

if any(infl)
    Wr = W; Wr(infl, :) = NaN;
    I.afterRemoval = rise_test_20r(Wr, xmin, sprintf('%s, influential flies removed', label));
    carried = sum(I.afterRemoval.RiseRatio <= I.afterRemoval.nullMedian) >= 3;
    if carried
        I.verdict = 'the rise is carried by a few flies';
    else
        I.verdict = 'some flies matter, but the rise does not depend on them alone';
    end
else
    I.verdict = 'no influential fly: the rise is not carried by any single fly';
end
fprintf('  influence verdict: %s\n', I.verdict);
end
