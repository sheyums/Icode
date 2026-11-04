function R = sleepVsCycleTests(T, options)
%SLEEPVSCYCLETESTS  Within-genotype Spearman of sleep against cycle direction and flux, Holm over 14 (PRESPEC 22; sqw22.py, sqw22h.py).
%
%   R = sleepVsCycleTests(T)                                      % observed level (PRESPEC 22)
%   R = sleepVsCycleTests(T, Direction="d", Seed=221, Secondary=false)   % hidden level (post hoc, sqw22h.py)
%
%   T   table, one row per fly x phase (SLEEPVSCYCLEROWS, or the Python's rows): genotype, fly,
%       phase (1-based), sleep (min/h), the direction column (c_obs or d), J, cycle (|J| >= 1e-12)
%   Direction  the direction outcome: "c_obs" (observed, default) or "d" (hidden)
%   NPerm      permutations per test (default 10000); NPermStrat for the stratified test (2000)
%   Seed       RandStream seed (default 22, as sqw22.py; 221 for sqw22h.py)
%   Secondary  also test CS DD and per0 LD light / dark (default true; outside the Holm family)
%
%   Per block: Spearman rho of sleep against the direction (rows with a finite value) and against J
%   (rows with a cycle only), permutation p (SPEARMANPERMTEST); Holm over the 14 PRIMARY tests (CS L, D;
%   per0 DD; fumin L, D; sss L, D x 2 outcomes); stratified Fisher z of the direction over the 7
%   primary blocks (STRATIFIEDSPEARMANTEST). rho is exact; p, Holm and the stratified p are Monte Carlo.
%   R   struct: blocks (table: genotype, phase, primary, n, sleep_med, n_dir, rho_dir, p_dir, holm_dir,
%       n_J, rho_J, p_J, holm_J, n_nocycle, sleep_nocycle, sleep_cycle), strat (z, rhoEquiv, p)
%
%   Example (from the Icode\sqw folder):
%       D = buildPerFlyData("D:\flies", fullfile("..", "hmm_audit", "sqw20b_fit2_theta0.5.jsonl"), Cache="perfly_cache.mat");
%       R = sleepVsCycleTests(sleepVsCycleRows(D)); disp(R.blocks)
%
%   See also SLEEPVSCYCLEROWS, SPEARMANPERMTEST, STRATIFIEDSPEARMANTEST, HOLMADJUST.

arguments
    T table
    options.Direction (1,1) string {mustBeMember(options.Direction, ["c_obs" "d"])} = "c_obs"
    options.NPerm (1,1) double = 10000
    options.NPermStrat (1,1) double = 2000
    options.Seed (1,1) double = 22
    options.Secondary (1,1) logical = true
end
rs = RandStream('twister', 'Seed', options.Seed); dv = options.Direction;
prim = {"CS", 1; "CS", 2; "per0DD", 1; "fumin", 1; "fumin", 2; "sss", 1; "sss", 2};
sec = {"CSDD", 1; "per0LD", 1; "per0LD", 2};
bl = prim; if options.Secondary, bl = [prim; sec]; end
B = table(); Xs = {}; Ys = {};
for k = 1:size(bl, 1)
    g = bl{k, 1}; p = bl{k, 2}; U = T(T.genotype == g & T.phase == p, :);
    if isempty(U), continue, end
    a = U(isfinite(U.(dv)), :); c = U(U.cycle, :); nc = U(~U.cycle, :);
    [rd, pd] = spearmanPermTest(a.sleep, a.(dv), options.NPerm, rs);
    [rj, pj] = spearmanPermTest(c.sleep, c.J, options.NPerm, rs);
    B = [B; table(g, p, k <= size(prim, 1), height(U), median(U.sleep), height(a), rd, pd, NaN, height(c), rj, pj, NaN, ...
        height(nc), median(nc.sleep), median(c.sleep), 'VariableNames', {'genotype', 'phase', 'primary', 'n', 'sleep_med', ...
        'n_dir', 'rho_dir', 'p_dir', 'holm_dir', 'n_J', 'rho_J', 'p_J', 'holm_J', 'n_nocycle', 'sleep_nocycle', 'sleep_cycle'})]; %#ok<AGROW>
    if k <= size(prim, 1), Xs{end + 1} = a.sleep; Ys{end + 1} = a.(dv); end %#ok<AGROW>
end
pr = find(B.primary); h = holmAdjust([B.p_dir(pr); B.p_J(pr)]);      % same family as sqw22: 7 x (dir, J)
B.holm_dir(pr) = h(1:numel(pr)); B.holm_J(pr) = h(numel(pr) + 1:end);
[z, pz, re] = stratifiedSpearmanTest(Xs, Ys, options.NPermStrat, rs);
R.blocks = B; R.strat = struct('z', z, 'rhoEquiv', re, 'p', pz); R.direction = dv;
end
