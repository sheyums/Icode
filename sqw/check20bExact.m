function R = check20bExact(dataRoot, jsonlFile, fitsFile, options)
%CHECK20BEXACT  Exact check of the port's HMM likelihood against cClaude's stored 20b fits (PORT_VALIDATION stage 2, Q2).
%
%   R = check20bExact(dataRoot, jsonlFile, fitsFile)
%   R = check20bExact(..., Tol=1e-6, EMFromStored=true)
%
%   dataRoot   folder with CS, per0LD, per0DD, fumin, sss (all 57 sss files: the
%              stored fits index sss in that set, A2 of ANSWERS_for_lClaude)
%   jsonlFile  hmm_audit/sqw20b_fit2_theta0.5.jsonl (one fly per line; fly is the
%              0-based index in sorted-glob file order)
%   fitsFile   sqw/validation/subset20b_fits.mat (the port's own 20-start fits,
%              from CHECKSUBSET20B); only these flies are checked
%
%   For each fly, three things:
%   1. Identity: T, nchg and the Q/W/S fractions of the port's symbols equal the
%      stored ones (the record is the same fly, read the same way).
%   2. Exact: HMMFORWARD at the stored (pi, P, B) for K = 2, 3, 4, and at the
%      labelled K = 3 fit (pi3, P3, B3), reproduces the stored ll to Tol.
%   3. Statistical: the port's own 20-start K = 3 optimum against the stored one
%      (different random starts; a higher port ll is legitimate). With
%      EMFromStored, Baum-Welch started AT the stored fit must not lower ll.
%
%   R  table, one row per fly.

arguments
    dataRoot (1,1) string {mustBeFolder}
    jsonlFile (1,1) string {mustBeFile}
    fitsFile (1,1) string {mustBeFile}
    options.Tol (1,1) double {mustBePositive} = 1e-6
    options.EMFromStored (1,1) logical = true
end
L = readlines(jsonlFile); L = L(strlength(L) > 0);
rec = cellfun(@jsondecode, cellstr(L));
recG = string({rec.g})'; recFly = [rec.fly]';
S = load(fitsFile, 'F', 'jobs'); F = S.F; jobs = S.jobs;

R = table();
for j = 1:numel(jobs)
    g = jobs(j).g;
    d = dir(fullfile(dataRoot, g, '**', '*.txt'));
    files = sort(string(fullfile({d.folder}, {d.name})));
    [~, nm, ext] = fileparts(files); nm = nm + ext;
    i0 = find(nm == jobs(j).file) - 1;                  % 0-based, as stored
    if numel(i0) ~= 1, error('sqw:file', 'File %s not found once under %s.', jobs(j).file, g); end
    r0 = rec(recG == g & recFly == i0);
    if numel(r0) ~= 1, error('sqw:record', 'No unique stored record for %s fly %d.', g, i0); end

    x = readCounts(files(i0 + 1));
    y = sqwSymbols(x, LightDark=jobs(j).ld);
    if jobs(j).ld, ph = minutePhases(numel(y)); nP = 2;
    else, ph = minutePhases(numel(y), Scheme="none"); nP = 1; end

    % 1. identity
    n = nnz(y >= 1 & y <= 3);
    frac = [nnz(y == 1), nnz(y == 2), nnz(y == 3)] / n;
    fracStored = [r0.frac.Q, r0.frac.W, r0.frac.S];
    nchg = nnz(diff(y) ~= 0);

    % 2. exact: forward at the stored parameters
    dK = nan(1, 3);
    for K = 2:4
        s = r0.F.(sprintf('x%d', K));
        [pi0, P, B] = storedParams(s.pi, s.P, s.B, K, nP);
        [~, ~, ll] = hmmForward(pi0, P, B, y, ph);
        dK(K - 1) = ll - s.ll;
    end
    [pi3, P3, B3] = storedParams(r0.pi3, r0.P3, r0.B3, 3, nP);
    [~, ~, ll3l] = hmmForward(pi3, P3, B3, y, ph);

    % 3. statistical: own optimum, and EM from the stored fit
    llStored = r0.F.x3.ll;
    gainEM = NaN; itEM = NaN;
    if options.EMFromStored
        [pi0, P, B] = storedParams(r0.F.x3.pi, r0.F.x3.P, r0.F.x3.B, 3, nP);
        st = struct('pi', pi0, 'P', P, 'B', B);
        fitS = hmmFitPhased(y, ph, nP, 3, Start=st);
        gainEM = fitS.ll - llStored; itEM = fitS.iterations;
    end

    r = table();
    r.genotype = g; r.file = jobs(j).file; r.fly0 = i0;
    r.T_port = numel(y); r.T_stored = r0.T;
    r.nchg_port = nchg; r.nchg_stored = r0.nchg;
    r.frac_maxdiff = max(abs(frac - fracStored));
    r.ll_stored_K3 = llStored;
    r.dll_K2 = dK(1); r.dll_K3 = dK(2); r.dll_K4 = dK(3);
    r.dll_K3_labelled = ll3l - llStored;
    r.exact_pass = all(abs([dK, ll3l - llStored]) <= options.Tol) && numel(y) == r0.T ...
        && nchg == r0.nchg && r.frac_maxdiff <= 1e-12;
    r.ll_port_own = F(j).ll; r.own_minus_stored = F(j).ll - llStored;
    r.stored_iterations = r0.F.x3.it; r.port_iterations = F(j).iterations;
    r.EM_from_stored_gain = gainEM; r.EM_from_stored_iterations = itEM;
    R = [R; r]; %#ok<AGROW>
end

fprintf('\nExact check (forward at stored parameters, Tol %.0e) and own optimum vs stored, K = 3\n', options.Tol);
for j = 1:height(R)
    fprintf('  %-6s fly %2d  T %d/%d  max|dll| K2-4 %.1e  labelled %.1e  %s   own-stored %+9.4f  EM-from-stored gain %+.2e\n', ...
        R.genotype(j), R.fly0(j), R.T_port(j), R.T_stored(j), max(abs([R.dll_K2(j) R.dll_K3(j) R.dll_K4(j)])), ...
        abs(R.dll_K3_labelled(j)), passWord(R.exact_pass(j)), R.own_minus_stored(j), R.EM_from_stored_gain(j));
end
fprintf('Exact: %d / %d flies pass. Own optimum within 1e-3 of stored: %d; above by more: %d; below by more: %d.\n', ...
    nnz(R.exact_pass), height(R), nnz(abs(R.own_minus_stored) <= 1e-3), ...
    nnz(R.own_minus_stored > 1e-3), nnz(R.own_minus_stored < -1e-3));
end

function [pi0, P, B] = storedParams(piS, PS, BS, K, nP)
% jsondecode gives P as nP x K x K (or K x K when nP = 1); 0-based state order is kept
pi0 = double(piS(:));
if nP == 1 && ismatrix(PS) && isequal(size(PS), [K K]), P = double(PS);
else, P = permute(double(reshape(PS, nP, K, K)), [2 3 1]); end
B = double(reshape(BS, K, 3));
end

function w = passWord(tf)
if tf, w = 'PASS'; else, w = 'FAIL'; end
end
