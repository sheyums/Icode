function D = planAData(dataRoot, jsonlFile, options)
%PLANADATA  Per fly x phase data behind Figs A1-A3 (PRESPEC 20s-A; sqw17.py), six genotypes in the paper's order.
%
%   D = planAData(dataRoot, jsonlFile)
%   D = planAData(dataRoot, jsonlFile, NPaths=200, NBoot=2000, Cache="")
%
%   dataRoot   one subfolder of count files per genotype: CS, CSDD, per0LD, per0DD,
%              fumin, sss54 (the 54 files of the 20s-A run) and sss (all 57; the
%              20b fits are indexed in that sorted list)
%   jsonlFile  hmm_audit/sqw20b_fit2_theta0.5.jsonl (READSTORED20BFITS)
%   Cache      .mat file: loaded if it exists, written otherwise
%
%   D  struct:
%     order     ["CS" "CSDD" "per0LD" "per0DD" "fumin" "sss"] (rule: CS -> CS DD ->
%               per0 LD -> per0 DD -> fumin -> sss; light before dark)
%     G         (1 x 6) SQWOBSERVABLEGROUP outputs in that order; perFly gains
%               hmm_deg, J_hmm, EPR_hmm, J_hmm_model, EPR_hmm_model (PLANAFLYHMM;
%               NaN and hmm_deg = true for CS DD, which has no stored fit)
%     hmmMatch  per genotype, how many flies were matched to a stored fit BY FILE NAME
%   The sss observable set is sss54 (as the 20s-A ledger); its HMM fits are found by
%   name in the 57-file list, never by position.
%
%   See also PLANAFLYHMM, EXPORTFIGA1, EXPORTFIGA2, EXPORTFIGA3.

arguments
    dataRoot (1,1) string {mustBeFolder}
    jsonlFile (1,1) string {mustBeFile}
    options.NPaths (1,1) double {mustBePositive, mustBeInteger} = 200
    options.NBoot (1,1) double {mustBePositive, mustBeInteger} = 2000
    options.Cache (1,1) string = ""
end
if strlength(options.Cache) > 0 && isfile(options.Cache)
    S = load(options.Cache, 'D'); D = S.D; fprintf('planAData: loaded %s\n', options.Cache); return
end
D.order = ["CS" "CSDD" "per0LD" "per0DD" "fumin" "sss"];
obsFolder = ["CS" "CSDD" "per0LD" "per0DD" "fumin" "sss54"];
fitFolder = ["CS" "" "per0LD" "per0DD" "fumin" "sss"];
dd = [false true false true false false];
fits = readStored20bFits(jsonlFile);
hmmCols = ["hmm_deg" "J_hmm" "EPR_hmm" "J_hmm_model" "EPR_hmm_model"];
for k = 1:6
    t0 = tic;
    G = sqwObservableGroup(fullfile(dataRoot, obsFolder(k)), D.order(k), LightDark=~dd(k), NBoot=options.NBoot);
    T = G.perFly; T.hmm_deg = true(height(T), 1);
    for c = hmmCols(2:end), T.(c) = NaN(height(T), 1); end
    nMatch = 0;
    if strlength(fitFolder(k)) > 0
        allF = sqwListFiles(fullfile(dataRoot, fitFolder(k)));
        [~, an, ae] = fileparts(allF); an = an + ae;
        Fg = fits([fits.g] == D.order(k));
        assert(numel(Fg) == numel(allF), 'planAData:fits', '%s: %d stored fits, %d files in %s.', ...
            D.order(k), numel(Fg), numel(allF), fitFolder(k));
        files = G.files;
        H = cell(numel(files), 1); hit = false(numel(files), 1);
        parfor i = 1:numel(files)
            [~, n, e] = fileparts(files(i)); j = find(an == n + e);
            if isscalar(j)
                f = Fg([Fg.fly0] == j - 1);
                [y, ph, nP] = sqwFlyInput(allF(j), ~dd(k));
                assert(numel(y) == f.T, 'planAData:T', '%s: %d minutes, stored fit has %d.', n + e, numel(y), f.T);
                H{i} = planAFlyHMM(y, ph, nP, f, NPaths=options.NPaths, RandomSeed=20 + 1000 * k + (j - 1));
                hit(i) = true;
            end
        end
        nMatch = nnz(hit);
        for i = find(hit)'
            r = find(T.fly == i);
            for c = hmmCols, T.(c)(r) = [H{i}.(c)]'; end
        end
    end
    G.perFly = T;
    D.G(k) = G; D.hmmMatch.(D.order(k)) = nMatch;
    fprintf('planAData: %-6s %2d flies (%d matched to a stored fit by name), %.0f s\n', D.order(k), ...
        numel(G.files), nMatch, toc(t0));
end
if strlength(options.Cache) > 0, save(options.Cache, 'D'); end
end
