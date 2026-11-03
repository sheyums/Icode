function R = checkHmmRefitSubset(dataRoot, csvFile, outFile, options)
%CHECKHMMREFITSUBSET  Refit a subset of flies (20b, K = 3) and compare hidden round trips with cClaude's per-fly table.
%
%   R = checkHmmRefitSubset(dataRoot, csvFile, outFile)
%   R = checkHmmRefitSubset(..., FliesPerGenotype=2, NPaths=200, Starts=20)
%
%   dataRoot  folder with CS, per0LD, per0DD, fumin, sss54 (the files cClaude used)
%   csvFile   hmm_audit/roundtrip_per_fly_theta0.5.csv (20g, from cClaude's 20b fits)
%   outFile   .mat for the fits (kept: the exact log-likelihood check needs them
%             when cClaude's stored fits arrive)
%
%   Statistical, not exact: for each fly x phase x level (hidden_all, hidden_matched)
%   z = (port - Python) / (sqrt(2) * SD over paths / sqrt(NPaths)) for p_quiet_entry,
%   p_quiet_exit and the four fractions. If the port reached the same optimum, the
%   z values behave like N(0,1). A different optimum shows up as large |z| on a fly.
%   Flies run in parallel (parfor) when a pool is available.

arguments
    dataRoot (1,1) string {mustBeFolder}
    csvFile (1,1) string {mustBeFile}
    outFile (1,1) string
    options.FliesPerGenotype (1,1) double {mustBePositive, mustBeInteger} = 2
    options.NPaths (1,1) double {mustBePositive, mustBeInteger} = 200
    options.Starts (1,1) double {mustBePositive, mustBeInteger} = 20
end
T = readtable(csvFile, 'TextType', 'string');
geno = ["CS" "per0LD" "per0DD" "fumin" "sss"]; folder = ["CS" "per0LD" "per0DD" "fumin" "sss54"];
jobs = struct('g', {}, 'file', {}, 'path', {}, 'ld', {});
for k = 1:numel(geno)
    d = dir(fullfile(dataRoot, folder(k), '**', '*.txt'));
    f = sort(string(fullfile({d.folder}, {d.name})));
    [~, nm, ext] = fileparts(f); nm = nm + ext;
    for i = 1:options.FliesPerGenotype
        jobs(end + 1) = struct('g', geno(k), 'file', nm(i), 'path', f(i), 'ld', geno(k) ~= "per0DD"); %#ok<AGROW>
    end
end
n = numel(jobs); out = cell(n, 1); nPaths = options.NPaths; nStarts = options.Starts;
parfor j = 1:n
    out{j} = oneFly(jobs(j), nPaths, nStarts);
end
F = vertcat(out{:});
save(outFile, 'F', 'jobs', '-v7');

stats = ["p_quiet_entry" "p_quiet_exit" "frac_forward" "frac_reverse" "frac_direct" "frac_quiet_both"];
R = table();
for j = 1:n
    f = F(j);
    for p = 1:numel(f.phaseNames)
        for lev = ["hidden_all" "hidden_matched"]
            row = T(T.genotype == jobs(j).g & T.file == jobs(j).file & T.phase == f.phaseNames(p) & T.level == lev, :);
            if height(row) ~= 1, continue, end
            S = f.(lev){p};
            if isempty(S), continue, end
            for s = stats
                r = table();
                r.genotype = jobs(j).g; r.file = jobs(j).file; r.phase = f.phaseNames(p); r.level = lev;
                r.stat = s; r.port = S.(s); r.python = row.(s); r.se = S.(s + "_se");
                % binomial floor: when every path agrees the path SD is 0, but one path in
                % NPaths differing is still sampling noise (e.g. 1 vs 0.995 at 200 paths)
                pp = min(max((r.port + r.python) / 2, 0), 1);
                if r.se == 0, r.se = sqrt(pp * (1 - pp) / nPaths); end   % only where all paths agree
                r.z = (r.port - r.python) / (sqrt(2) * max(r.se, 1e-6));
                R = [R; r]; %#ok<AGROW>
            end
        end
    end
end
fprintf('\nPer fly: log-likelihood (K = 3), degenerate H1, max |z| over 6 statistics x phases x levels\n');
for j = 1:n
    zz = abs(R.z(R.file == jobs(j).file));
    fprintf('  %-6s %-60s ll %.4f  deg %d  max|z| %.2f\n', jobs(j).g, jobs(j).file, F(j).ll, F(j).deg, max(zz));
end
fprintf('All comparisons: n %d, |z| > 2 in %.3f (expect ~0.05), |z| > 3 in %.3f, max |z| %.2f\n', ...
    height(R), mean(abs(R.z) > 2), mean(abs(R.z) > 3), max(abs(R.z)));
end

function f = oneFly(job, nPaths, nStarts)
x = readCounts(job.path);
y = sqwSymbols(x, LightDark=job.ld);
if job.ld, ph = minutePhases(numel(y)); nP = 2; f.phaseNames = ["light"; "dark"];
else, ph = minutePhases(numel(y), Scheme="none"); nP = 1; f.phaseNames = "whole"; end
fit = hmmFitPhased(y, ph, nP, 3, Starts=nStarts);
[P, B, pi0] = labelStates(fit.P, fit.B, fit.pi);
f.ll = fit.ll; f.P = P; f.B = B; f.pi = pi0; f.deg = isDegenerateH1(B); f.iterations = fit.iterations;
Z = hmmFFBS(pi0, P, B, y, ph, NPaths=nPaths, RandomSeed=1);
hours = accumarray(ph, 1, [nP 1]) / 60;
f.hidden_all = cell(nP, 1); f.hidden_matched = cell(nP, 1);
for lev = ["hidden_all" "hidden_matched"]
    H = zeros(nPaths, 4, nP);
    for k = 1:nPaths
        if lev == "hidden_all", c = roundTripsHidden(Z(k, :)', ph, nP);
        else, c = roundTripsHidden(Z(k, :)', ph, nP, MinSleep=5, MinGap=10); end
        H(k, :, :) = reshape(c', 1, 4, nP);
    end
    for p = 1:nP
        C = H(:, :, p); s = roundTripStats(C, hours(p));
        if isempty(s), f.(lev){p} = []; continue, end      % no round trip on any path
        nn = sum(C, 2); ok = nn > 0; Fr = C(ok, :) ./ nn(ok);
        qi = Fr(:, 1) + Fr(:, 4); qo = Fr(:, 2) + Fr(:, 4);
        sd = [std(qi), std(qo), std(Fr, 0, 1)] / sqrt(nnz(ok));
        S = struct('p_quiet_entry', s.p_q_entry, 'p_quiet_exit', s.p_q_exit, 'frac_forward', s.frac(1), ...
            'frac_reverse', s.frac(2), 'frac_direct', s.frac(3), 'frac_quiet_both', s.frac(4));
        nm = fieldnames(S);
        for q = 1:numel(nm), S.([nm{q} '_se']) = sd(q); end
        f.(lev){p} = S;
    end
end
end
