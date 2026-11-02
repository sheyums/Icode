% MATLAB counterpart of sqw_r3.py: rule-of-three medians (ledger amendment 2026-11-02)
assert(isequal(fileparts(which('modelEnergeticsR3')), pwd), 'run from Icode/sqw');
sp = fileparts(mfilename('fullpath')); C = load(fullfile(sp, 'out_figA', 'planAData_cache.mat')); f = fieldnames(C); D = C.(f{1});
fits = readStored20bFits(fullfile(pwd, '..', 'hmm_audit', 'sqw20b_fit2_theta0.5.jsonl'));
fitFolder = ["CS" "" "per0LD" "per0DD" "fumin" "sss"]; rows = {};
c = parcluster('Processes'); c.NumWorkers = 8; delete(gcp('nocreate')); parpool(c, 8);
for k = 1:6
    G = D.G(k); T = G.perFly; g = D.order(k);
    if strlength(fitFolder(k)) > 0
        allF = sqwListFiles(fullfile(sp, 'data_cc', fitFolder(k))); [~, an, ae] = fileparts(allF); an = an + ae; Fg = fits([fits.g] == g);
    end
    Hn = cell(numel(G.flies), 1); fl = G.flies;
    parfor i = 1:numel(fl)                    % mean posterior-path pair counts, labelled K = 3 stored fit
        if strlength(fitFolder(k)) == 0, continue, end
        [~, n, e] = fileparts(fl(i).file); j = find(an == n + e); if ~isscalar(j), continue, end
        F = Fg([Fg.fly0] == j - 1); if F.deg, continue, end
        Z = hmmFFBS(F.pi(:), F.P, F.B, fl(i).y, fl(i).ph, NPaths=200, RandomSeed=7000 + 1000 * k + j);
        Ns = zeros(3, 3, fl(i).nP); for z = 1:size(Z, 1), Ns = Ns + pairCounts(Z(z, :)', fl(i).ph, fl(i).nP); end
        Hn{i} = struct('P', F.P, 'N', Ns / size(Z, 1));
    end
    for r = 1:height(T)
        if ~T.refOK(r), continue, end
        i = T.fly(r); p = T.phase(r); P = fl(i).chain.P(:, :, p); N = reshape(T.N(r, :), 3, 3);
        e0 = modelEnergetics(P); e = modelEnergeticsR3(P, N); L = energyLandscape(e.Pr3);
        row = {g, i, p, T.hasCycle(r), e0.A, e0.EPR, e.A, e.EPR, e.r3, L.E2_minus_E1, L.E3_minus_E1, NaN, NaN, false};
        if ~isempty(Hn{i})
            eh = modelEnergeticsR3(Hn{i}.P(:, :, p), Hn{i}.N(:, :, p));
            row(12:14) = {T.EPR_hmm_model(r), eh.EPR, eh.r3};
        end
        rows(end + 1, :) = row; %#ok<AGROW>
    end
    fprintf('%s done\n', g);
end
R = cell2table(rows, 'VariableNames', {'g', 'fly', 'p', 'cycle', 'A_old', 'EPR_old', 'A_r3', 'EPR_r3', 'r3', 'E21_r3', 'E31_r3', 'EPRh_old', 'EPRh_r3', 'r3h'});
writetable(R, fullfile(sp, 'r3_port.csv'));
fprintf('genotype phase | n cycle | r3 | A old -> r3 | chain EPR old -> r3 | HMM model EPR old -> r3 (r3 used)\n');
for g = D.order
    for p = unique(R.p(R.g == g))'
        rr = R(R.g == g & R.p == p & R.cycle, :); hh = R(R.g == g & R.p == p & ~isnan(R.EPRh_r3), :);
        fprintf('%-6s %d | %2d | %2d | %+.2f -> %+.2f | %.3f -> %.3f | %.3f -> %.3f (%d)\n', g, p - 1, height(rr), nnz(rr.r3), ...
            pm(rr.A_old), pm(rr.A_r3), pm(rr.EPR_old), pm(rr.EPR_r3), pm(hh.EPRh_old), pm(hh.EPRh_r3), nnz(hh.r3h));
    end
end
fprintf('observed fly-phases with a cycle still non-finite after r3: %d; observed r3 used %d, hidden %d\n', ...
    nnz(R.cycle & ~(isfinite(R.A_r3) & isfinite(R.EPR_r3))), nnz(R.r3), nnz(R.r3h));
function m = pm(v)
% median as sqw_r3.py: non-finite (not NaN) ranked as +Inf, NaN dropped
v = v(~isnan(v)); v(~isfinite(v)) = Inf; if isempty(v), m = NaN; else, m = median(v); end
end
