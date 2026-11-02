function H = planARouteHidden(dataRoot, jsonlFile, options)
%PLANAROUTEHIDDEN  Hidden-level 20k inputs for R1 (sqw10.py): per-fly hidden D, its simulation null, and path records.
%
%   H = planARouteHidden(dataRoot, jsonlFile)
%   H = planARouteHidden(dataRoot, jsonlFile, NPaths=200, NSim=400, NPathSim=20, NRec=20, Cache="", Genotypes=...)
%
%   dataRoot   count files per genotype (CS, fumin, per0LD, per0DD, sss = all 57); jsonlFile the stored
%              20b fits. Flies with a NON-degenerate fit, as sqw10.py.
%   Per fly, with its labelled K = 3 fit (empty rows filled for simulation, as everywhere in the port):
%     hidD      (1 x nP) mean over NPaths FFBS paths of ROUTEDURATIOND(ROUTEDURATIONHIDDEN(z), R1, p),
%               NaN unless at least half the paths give a finite D
%     sims      (NSim x nP) the same statistic on NSim records simulated from the fit (NOT masked,
%               as sqw10.py) and decoded with NPathSim FFBS paths each
%     records   the ROUTEDURATIONHIDDEN records of the first NRec real paths, pooled (panel c)
%     Fdark     (F_SW + 1)/2 of the fit's last-phase matrix if it has a cycle (sqw11.py), else NaN
%   H  table, one row per fly. FFBS draws come from MATLAB streams: compare statistically.
%
%   See also ROUTEDURATIONHIDDEN, HMMFFBS, EXPORTFIGC7SI.

arguments
    dataRoot (1,1) string {mustBeFolder}
    jsonlFile (1,1) string {mustBeFile}
    options.NPaths (1,1) double {mustBePositive, mustBeInteger} = 200
    options.NSim (1,1) double {mustBeNonnegative, mustBeInteger} = 400
    options.NPathSim (1,1) double {mustBePositive, mustBeInteger} = 20
    options.NRec (1,1) double {mustBePositive, mustBeInteger} = 20
    options.Cache (1,1) string = ""
    options.Genotypes (1,:) string = ["CS" "fumin" "per0LD" "per0DD" "sss"]
    options.MaxFlies (1,1) double = Inf
end
if strlength(options.Cache) > 0 && isfile(options.Cache)
    C = load(options.Cache, 'H'); H = C.H; fprintf('planARouteHidden: loaded %s\n', options.Cache); return
end
fits = readStored20bFits(jsonlFile); H = table(); o = options;
for gk = 1:numel(o.Genotypes)
    g = o.Genotypes(gk); files = sqwListFiles(fullfile(dataRoot, g)); Fg = fits([fits.g] == g);
    use = find(~arrayfun(@(i) Fg([Fg.fly0] == i - 1).deg, (1:numel(files))'));
    use = use(1:min(numel(use), o.MaxFlies)); ld = g ~= "per0DD"; rows = cell(numel(use), 1); t0 = tic;
    parfor j = 1:numel(use)
        i = use(j); f = Fg([Fg.fly0] == i - 1); f.P = fillEmptyRows(f.P);
        [y, ph, nP, refs, names] = sqwFlyInput(files(i), ld);
        Z = hmmFFBS(f.pi(:), f.P, f.B, y, ph, NPaths=o.NPaths, RandomSeed=1000 * gk + i);
        [hd, rec] = hidD(Z, ph, nP, o.NRec);
        rs = RandStream('twister', 'Seed', 500000 + 1000 * gk + i); S = NaN(o.NSim, nP);
        for s = 1:o.NSim
            [~, ys] = hmmSimulatePhased(f.pi(:), f.P, f.B, ph, rs);
            Zs = hmmFFBS(f.pi(:), f.P, f.B, ys, ph, NPaths=o.NPathSim, RandomSeed=randi(rs, 2^30));
            S(s, :) = hidD(Zs, ph, nP, 0);
        end
        Fd = NaN; Pd = f.P(:, :, end);
        if isfinite(refs.(names(end))) && abs(modelEnergetics(Pd).J) >= 1e-12
            L = energyLandscape(Pd); Fd = (L.F_SW + 1) / 2;
        end
        [~, nm, ex] = fileparts(files(i));
        rows{j} = table(g, i - 1, string(nm) + ex, nP, {hd}, {S}, {rec}, Fd, ...
            'VariableNames', {'genotype', 'fly0', 'file', 'nP', 'hidD', 'sims', 'records', 'Fdark'});
    end
    H = [H; vertcat(rows{:})]; %#ok<AGROW>
    fprintf('planARouteHidden: %-6s %d flies, %.0f s\n', g, numel(use), toc(t0));
end
if strlength(options.Cache) > 0, save(options.Cache, 'H', '-v7.3'); end
end

function [d, rec] = hidD(Z, ph, nP, nRec)
n = size(Z, 1); v = NaN(n, nP); rec = [];
for k = 1:n
    R = routeDurationHidden(Z(k, :)', ph);
    for p = 1:nP, v(k, p) = routeDurationD(R, "R1_exit_vs_sleep", p); end
    if k <= nRec, rec = [rec; R]; end %#ok<AGROW>
end
d = NaN(1, nP);
for p = 1:nP
    if mean(isfinite(v(:, p))) >= 0.5, d(p) = mean(v(isfinite(v(:, p)), p)); end
end
end
