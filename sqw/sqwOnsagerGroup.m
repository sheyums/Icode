function O = sqwOnsagerGroup(G, options)
%SQWONSAGERGROUP  EPR against J and the Onsager conductance L per genotype (PRESPEC 20s-A, A2; sqw17_report).
%
%   O = sqwOnsagerGroup(G)
%   O = sqwOnsagerGroup(G, Points="count", Reference="CS", NBoot=2000, RandomSeed=20, OutDir="")
%
%   G          array of SQWOBSERVABLEGROUP outputs, one per genotype
%   Points     "count": count J and count EPR (rule of three);
%              "chain": the observed chain's exact J and EPR
%   Reference  genotype the others are compared with
%   WaldSet    genotypes entering the Wald test and its BH family (default:
%              all but Reference; the ledger used per0LD, per0DD, fumin, sss)
%
%   Points are fly x phase (both phases pooled), only fly-phases with a Q/W
%   reference, a chain with a cycle, and finite J and EPR.
%
%   O  table, one row per genotype:
%      n, nFly, a / rss1 (M1: EPR = b + a|J|), c / rss2 (M2: EPR = b + cJ^2),
%      better ("M1 linear" / "M2 quadratic"), shareM1 (bootstrap share where M1
%      wins), exponent and its bootstrap interval, L = 1/s through the origin,
%      se_s, L_lo / L_hi (fly bootstrap, type-7 percentiles), L_t_lo / L_t_hi
%      (ONSAGERFIT t interval), and against Reference: Wald z and p on 1/L,
%      BH-adjusted p, and the fly-bootstrap difference interval and p
%      (BOOTSTRAPDIFFERENCEP, the Methods' procedure).
%   O.Properties.UserData.points  the points; .boot  the bootstrap L vectors.
%   Wald is known to overstate (20j note); the bootstrap difference is the
%   comparison the Methods names.
%
%   See also ONSAGERFIT, EPRJMODELS, ONSAGERBOOTSTRAP, PLOTSQWONSAGER.

arguments
    G struct
    options.Points (1,1) string {mustBeMember(options.Points, ["count" "chain"])} = "count"
    options.Reference (1,1) string = "CS"
    options.WaldSet string = string.empty
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 20
    options.OutDir (1,1) string = ""
end
if options.Points == "count", kJ = "J"; kE = "EPR"; else, kJ = "J_chain"; kE = "EPR_chain"; end
rs = RandStream('twister', 'Seed', options.RandomSeed);
O = table(); pts = table(); boot = struct();
for g = 1:numel(G)
    T = G(g).perFly;
    T = T(T.refOK & T.hasCycle & isfinite(T.(kJ)) & isfinite(T.(kE)), :);
    J = T.(kJ); E = T.(kE); fly = T.fly;
    M = eprJModels(J, E);
    F = onsagerFit(J, E);
    [Lb, mb, nb] = bootAll(fly, J, E, options.NBoot, rs);
    r = table();
    r.genotype = G(g).genotype; r.n = numel(J); r.nFly = numel(unique(fly));
    r.a = M.a; r.rss1 = M.rss_linear; r.c = M.c; r.rss2 = M.rss_quadratic;
    if M.rss_linear < M.rss_quadratic, r.better = "M1 linear"; else, r.better = "M2 quadratic"; end
    r.shareM1 = mean(mb);
    r.exponent = M.exponent; e = nb(isfinite(nb));
    q = quantile7(e, [0.025 0.975]); r.exp_lo = q(1); r.exp_hi = q(2);
    r.L = F.L; r.s = F.s; r.se_s = F.se_s;
    q = quantile7(Lb, [0.025 0.975]); r.L_lo = q(1); r.L_hi = q(2);
    r.L_t_lo = F.ci_L(1); r.L_t_hi = F.ci_L(2);
    O = [O; r]; %#ok<AGROW>
    boot.(matlab.lang.makeValidName(G(g).genotype)) = Lb;
    T.genotype = repmat(G(g).genotype, height(T), 1);
    pts = [pts; T(:, ["genotype" "fly" "file" "phaseName" kJ kE])]; %#ok<AGROW>
end

% comparisons with the reference genotype
iR = find(O.genotype == options.Reference, 1);
O.wald_z = NaN(height(O), 1); O.wald_p = O.wald_z; O.wald_BH = O.wald_z;
O.diff_lo = O.wald_z; O.diff_hi = O.wald_z; O.diff_p = O.wald_z;
if ~isempty(iR)
    k = setdiff(1:height(O), iR);
    kw = k;
    if ~isempty(options.WaldSet), kw = find(ismember(O.genotype, options.WaldSet))'; kw = setdiff(kw, iR); end
    O.wald_z(kw) = (O.s(kw) - O.s(iR)) ./ hypot(O.se_s(kw), O.se_s(iR));
    O.wald_p(kw) = erfc(abs(O.wald_z(kw)) / sqrt(2));
    O.wald_BH(kw) = bh(O.wald_p(kw));
    LbR = boot.(matlab.lang.makeValidName(O.genotype(iR)));
    for j = k
        [ci, p] = bootstrapDifferenceP(boot.(matlab.lang.makeValidName(O.genotype(j))), LbR);
        O.diff_lo(j) = ci(1); O.diff_hi(j) = ci(2); O.diff_p(j) = p;
    end
end
O.Properties.UserData = struct('points', pts, 'boot', boot, 'pointType', options.Points);
if strlength(options.OutDir) > 0
    if ~isfolder(options.OutDir), mkdir(options.OutDir); end
    writetable(O, fullfile(options.OutDir, "onsager_" + options.Points + "_summary.csv"));
    writetable(pts, fullfile(options.OutDir, "onsager_" + options.Points + "_points.csv"));
end
end

function [Lb, m1, ex] = bootAll(fly, J, E, nBoot, rs)
% Flies resampled with all their points; L, model choice and exponent from each resample.
[~, ~, g] = unique(fly); nF = max(g);
rows = accumarray(g, (1:numel(g))', [nF 1], @(v) {v});
Lb = zeros(nBoot, 1); m1 = false(nBoot, 1); ex = NaN(nBoot, 1);
for b = 1:nBoot
    r = vertcat(rows{randi(rs, nF, nF, 1)});
    x = J(r) .^ 2; e = E(r);
    Lb(b) = (x' * x) / (x' * e);
    M = eprJModels(J(r), E(r));
    m1(b) = M.rss_linear < M.rss_quadratic; ex(b) = M.exponent;
end
end

function q = bh(p)
% Benjamini-Hochberg adjusted p-values (step-up, monotone).
[ps, o] = sort(p(:)); m = numel(ps);
a = min(1, cummin(ps .* m ./ (1:m)', 'reverse'));
q = zeros(m, 1); q(o) = a;
end
