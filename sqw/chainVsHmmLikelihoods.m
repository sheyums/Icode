function L = chainVsHmmLikelihoods(dataRoot, jsonlFile, options)
%CHAINVSHMMLIKELIHOODS  Per-fly log-likelihoods of PRESPEC 20n and 20o (sqw13.py, sqw14.py) for Fig A4 (a, b).
%
%   L = chainVsHmmLikelihoods(dataRoot, jsonlFile)
%   L = chainVsHmmLikelihoods(dataRoot, jsonlFile, Genotypes=["CS" "fumin"], Starts=20, Cache="")
%
%   dataRoot   one subfolder of count files per genotype (the 20b fits are indexed
%              in each folder's sorted file list, READSTOREDHMMFITS)
%   jsonlFile  hmm_audit/sqw20b_fit2_theta0.5.jsonl
%   Cache      .mat file: loaded if it exists, written otherwise
%
%   L  table, one row per fly: genotype, fly0, file, T (minutes), deg, and the full-record
%      log-likelihoods (missing minutes marginalised), 12-h light/dark or four 6-h phases:
%        MC1_12  first-order chain (CHAINFITPHASED)            k 14
%        HMM_12  the stored 20b K = 3 fit's ll (its F.x3.ll)    k 20
%        MC2_12  second-order chain (CHAIN2FITPHASED)          k 44
%        MC1_6   first-order chain, 6-h phases                 k 26
%        HMM_6   K = 3 HMM, 6-h phases, HMMFITPHASED 20 starts  k 32
%      dBIC_T = BIC(MC1_12) - BIC(HMM_12) with n = T (20n; > 0: BIC prefers the HMM),
%      gain = HMM_12 - MC1_12, R_mem = (MC2_12 - MC1_12) / gain,
%      R_drift = 1 - (HMM_6 - MC1_6) / gain  (20o; read only where gain > 3 nats)
%   HMM_6 comes from random starts drawn by MATLAB's stream, not numpy's: it can
%   land on a different (legitimately higher or lower) optimum than sqw14.py.
%
%   Example (from the Icode\sqw folder):
%       dataRoot = "D:\flies";          % one subfolder of 1-s count files per genotype
%       jsonl = fullfile("..", "hmm_audit", "sqw20b_fit2_theta0.5.jsonl");
%       L = chainVsHmmLikelihoods(dataRoot, jsonl, Genotypes=["CS" "fumin"], Cache="likelihoods.mat");
%       mean(L.dBIC_T > 0)   % share of flies where BIC prefers the HMM
%
%   See also EXPORTHIDDENVSOBSERVEDFIGURE, CHAINFITPHASED, CHAIN2FITPHASED, HMMFITPHASED.

arguments
    dataRoot (1,1) string {mustBeFolder}
    jsonlFile (1,1) string {mustBeFile}
    options.Genotypes (1,:) string = ["CS" "fumin"]
    options.Starts (1,1) double {mustBePositive, mustBeInteger} = 20
    options.Cache (1,1) string = ""
end
if strlength(options.Cache) > 0 && isfile(options.Cache)
    S = load(options.Cache, 'L'); L = S.L; fprintf('chainVsHmmLikelihoods: loaded %s\n', options.Cache); return
end
fits = readStoredHmmFits(jsonlFile); L = table(); nSt = options.Starts;
for k = 1:numel(options.Genotypes)
    g = options.Genotypes(k); files = sqwListFiles(fullfile(dataRoot, g));
    Fg = fits([fits.g] == g);
    assert(numel(Fg) == numel(files), 'chainVsHmmLikelihoods:fits', '%s: %d fits, %d files.', g, numel(Fg), numel(files));
    rows = cell(numel(files), 1); t0 = tic;
    parfor i = 1:numel(files)
        f = Fg([Fg.fly0] == i - 1);
        [y, ph] = sqwFlyInput(files(i), true);
        assert(numel(y) == f.T, 'chainVsHmmLikelihoods:T', '%s fly %d: %d minutes, fit %d.', g, i - 1, numel(y), f.T);
        ph6 = minutePhases(numel(y), Scheme="6h");
        r = struct();
        r.MC1_12 = chainFitPhased(y, ph, 2);
        r.HMM_12 = f.F.x3.ll;
        r.MC2_12 = chain2FitPhased(y, ph, 2);
        r.MC1_6 = chainFitPhased(y, ph6, 4);
        h6 = hmmFitPhased(y, ph6, 4, 3, Starts=nSt, RandomSeed=20 + 1000 * k + (i - 1));
        r.HMM_6 = h6.ll;
        [~, nm, ex] = fileparts(files(i));
        rows{i} = [table(g, i - 1, string(nm) + ex, numel(y), logical(f.deg), ...
            'VariableNames', {'genotype', 'fly0', 'file', 'T', 'deg'}), struct2table(r)];
    end
    L = [L; vertcat(rows{:})]; %#ok<AGROW>
    fprintf('chainVsHmmLikelihoods: %-6s %2d flies, %.0f s\n', g, numel(files), toc(t0));
end
lnT = log(L.T);
L.dBIC_T = (-2 * L.MC1_12 + 14 * lnT) - (-2 * L.HMM_12 + 20 * lnT);
L.gain = L.HMM_12 - L.MC1_12;
L.R_mem = (L.MC2_12 - L.MC1_12) ./ L.gain;
L.R_drift = 1 - (L.HMM_6 - L.MC1_6) ./ L.gain;
if strlength(options.Cache) > 0, save(options.Cache, 'L'); end
end
