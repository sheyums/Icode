function T = sleepVsCycleRows(D, options)
%SLEEPVSCYCLEROWS  Per fly x phase sleep, cycle direction and flux, observed or hidden (rows of sqw22.py / sqw22h.py).
%
%   T = sleepVsCycleRows(D)                                              % observed level
%   T = sleepVsCycleRows(D, Level="hidden", DataRoot=..., JsonlFile=..., NPaths=200, RandomSeed=221)
%
%   D   BUILDPERFLYDATA output
%   Observed (PRESPEC 22), every fly x phase with a Q/W reference:
%     sleep = 60 x the fraction of the phase's SCORED minutes that are S (min/h); c_obs (COBSMINUTES);
%     J = J_chain of the FITTED chain (not the rule-of-three matrix); cycle = |J_chain| >= 1e-12, else J = 0.
%   Hidden (post hoc, sqw22h.py), fly x phases with a reference and a non-degenerate stored 20b fit
%   (not CS DD, which was not fitted):
%     sleep = 60 x pi_S, the stationary probability of the sleep-like state of the fly's labelled K = 3
%     matrix for that phase; J = J_hmm (D: count J averaged over posterior paths); cycle = |J| >= 1e-12;
%     d = (N_SW - N_WS) / (N_QS + N_WS), net forward turns per sleep onset, computed on EACH of NPaths
%     FFBS paths and averaged over the paths with an onset, defined only if at least half have one.
%   DataRoot   folder with the count-file subfolders CS, per0LD, per0DD, fumin, sss (all 57: the stored
%              fits index sss in that set); needed for Level="hidden" to match fits to flies by name
%   JsonlFile  hmm_audit/sqw20b_fit2_theta0.5.jsonl
%   T   table: genotype, fly (index into D.G(k).flies), phase (1-based), sleep, c_obs (observed) or d
%       (hidden), J, cycle
%
%   Example (from the Icode\sqw folder):
%       jsonl = fullfile("..", "hmm_audit", "sqw20b_fit2_theta0.5.jsonl");
%       D = buildPerFlyData("D:\flies", jsonl, Cache="perfly_cache.mat");
%       Ro = sleepVsCycleTests(sleepVsCycleRows(D));
%       Rh = sleepVsCycleTests(sleepVsCycleRows(D, Level="hidden", DataRoot="D:\flies", JsonlFile=jsonl), ...
%                              Direction="d", Seed=221, Secondary=false);
%
%   See also SLEEPVSCYCLETESTS, BUILDPERFLYDATA, HMMFFBS, PAIRCOUNTS.

arguments
    D (1,1) struct
    options.Level (1,1) string {mustBeMember(options.Level, ["observed" "hidden"])} = "observed"
    options.DataRoot (1,1) string = ""
    options.JsonlFile (1,1) string = ""
    options.NPaths (1,1) double {mustBePositive, mustBeInteger} = 200
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 221
end
T = table(); hid = options.Level == "hidden"; fits = []; Fg = []; an = strings(0, 1);
if hid, fits = readStoredHmmFits(options.JsonlFile); end
fitFolder = containers.Map(["CS" "per0LD" "per0DD" "fumin" "sss"], ["CS" "per0LD" "per0DD" "fumin" "sss"]);
for k = 1:numel(D.G)
    g = D.order(k); G = D.G(k); P = G.perFly;
    if hid
        if ~isKey(fitFolder, g), continue, end
        allF = sqwListFiles(fullfile(options.DataRoot, fitFolder(g))); [~, an, ae] = fileparts(allF); an = an + ae;
        Fg = fits([fits.g] == g);
    end
    rows = cell(numel(G.flies), 1); fl = G.flies;
    parfor i = 1:numel(fl)
        f = fl(i); U = P(P.fly == i & P.refOK, :); r = table();
        if isempty(U), continue, end
        if ~hid
            for q = 1:height(U)
                p = U.phase(q); m = f.ph == p & f.y >= 1 & f.y <= 3; J = U.J_chain(q); cyc = isfinite(J) && abs(J) >= 1e-12; if ~cyc, J = 0; end
                r = [r; table(g, i, p, 60 * mean(f.y(m) == 3), U.c_obs(q), J, cyc, ...
                    'VariableNames', {'genotype', 'fly', 'phase', 'sleep', 'c_obs', 'J', 'cycle'})]; %#ok<AGROW>
            end
        else
            [~, n, e] = fileparts(f.file); j = find(an == n + e);
            if ~isscalar(j), continue, end
            F = Fg([Fg.fly0] == j - 1); if F.deg, continue, end
            Z = hmmFFBS(F.pi(:), F.P, F.B, f.y, f.ph, NPaths=options.NPaths, RandomSeed=options.RandomSeed + 1000 * k + j);
            Np = zeros(3, 3, f.nP, size(Z, 1));
            for z = 1:size(Z, 1), Np(:, :, :, z) = pairCounts(Z(z, :)', f.ph, f.nP); end
            for q = 1:height(U)
                p = U.phase(q); J = U.J_hmm(q); if ~isfinite(J), continue, end
                pis = stationary_dist(F.P(:, :, p)); on = squeeze(Np(1, 3, p, :) + Np(2, 3, p, :)); net = squeeze(Np(3, 2, p, :) - Np(2, 3, p, :));
                dv = NaN; if mean(on > 0) >= 0.5, dv = mean(net(on > 0) ./ on(on > 0)); end
                r = [r; table(g, i, p, 60 * pis(3), dv, J, abs(J) >= 1e-12, ...
                    'VariableNames', {'genotype', 'fly', 'phase', 'sleep', 'd', 'J', 'cycle'})]; %#ok<AGROW>
            end
        end
        rows{i} = r;
    end
    T = [T; vertcat(rows{:})]; %#ok<AGROW>
end
end
