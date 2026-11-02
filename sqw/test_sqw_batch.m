function test_sqw_batch()
%TEST_SQW_BATCH  Batch symbol and HMM tools (SQWBATCHSYMBOLS, SQWBATCHHMM, lab coding): stress tests on synthetic flies.
%
%   Run by name from the sqw folder:   cd <Icode>/sqw; test_sqw_batch
%
%    1. sqwToLabCode / labCodeToSqw: every code maps as documented (S 0, Q 1, W 2, missing NaN,
%       unknown wake -1) and back; anything else (3, 0.5, Inf, port 5) is refused
%    2. sqwBatchSymbols, LD: each written symbol column = sqwToLabCode(sqwSymbols(counts));
%       zt_h, phase, active_frac, sleep_s as SQWSYMBOLS' info; S+Q+W+missing = minutes overall
%       and per phase; r_light / r_dark = the references; outputs land BESIDE the input folder,
%       so a second run still finds exactly the count files
%    3. sqwBatchSymbols, DD: one phase, r_whole, no light / dark columns
%    4. file-list input without OutDir refused; a header line in a count file refused (READCOUNTS)
%    5. sqwBatchHMM reads the symbol files and the count files to the SAME log-likelihoods
%       (the two input routes are one analysis)
%    6. BIC columns = HMMBIC with n = minutes and n = symbol changes; bestK = argmin; flag =
%       some K ~= 3 below BIC(3) - 10; dBIC columns = differences
%    7. labelled K = 3: emission and transition rows sum to 1, H3 has the largest P(S), H2 the
%       larger P(W) of the other two, deg3 = ISDEGENERATEH1
%    8. dom_K2 / dom_K4: one symbol per state, sorted S, Q, W; nS + nQ + nW = 4
%    9. Resume: a cached fit is reused when Ks, Starts, seed and T match (a planted
%       log-likelihood comes back), and ignored when Starts differs
%   10. malformed symbol files refused: code 3, phase 3, missing phase column
%   11. recovery: 10 days simulated from a K = 3 model with distinct emissions -> bestK_T = 3,
%       and the fitted log-likelihood >= the log-likelihood at the true parameters
%   12. recovery: 10 days from a K = 2 model -> bestK_T = 2

here = fileparts(mfilename('fullpath'));
assert(isequal(fileparts(which('sqwBatchHMM')), here), 'Run from the sqw folder.');
nP = 0; nF = 0; tmp = tempname; mkdir(tmp); c0 = onCleanup(@() rmdir(tmp, 's'));
fprintf('\nRunning batch symbol / HMM tests...\n\n');
rs = RandStream('twister', 'Seed', 7);
fast = {'Starts', 4, 'MaxIter', 300, 'Plot', false};

%% 1. coding
try
    y = [1 2 3 0 4]'; s = sqwToLabCode(y);
    ok = isequaln(s, [1 2 0 NaN -1]') && isequal(labCodeToSqw(s), y);
    bad = 0;
    for v = {3, 0.5, Inf}, try, labCodeToSqw(v{1}); catch, bad = bad + 1; end, end
    try, sqwToLabCode(5); catch, bad = bad + 1; end
    [nP, nF] = rep(ok && bad == 4, nP, nF, 1, 'lab coding S 0 / Q 1 / W 2 / NaN / -1 and back; bad codes refused', ...
        sprintf('map ok %d, refused %d of 4', ok, bad));
catch e, [nP, nF] = rep(false, nP, nF, 1, 'lab coding', e.message); end

%% 2. sqwBatchSymbols LD
ld = fullfile(tmp, 'ld'); mkdir(ld);
for i = 1:2, writematrix(syntheticCounts(rs, 2), fullfile(ld, "fly" + i + ".txt")); end
try
    X = sqwBatchSymbols(ld); ok = true; why = "";
    for i = 1:2
        x = readCounts(fullfile(ld, "fly" + i + ".txt")); [y, refs, info] = sqwSymbols(x);
        t = readtable(fullfile(X.outDir, "fly" + i + "_symbols.csv")); T = numel(y);
        r = X.summary(i, :); m = (1:T)';
        ok = ok && isequaln(t.symbol, sqwToLabCode(y)) && isequal(t.minute, m) && max(abs(t.zt_h - mod((m - 1) / 60, 24))) < 1e-12 ...
            && isequal(t.phase, minutePhases(T)) && max(abs(t.active_frac - info.activeFraction(1:T))) < 1e-12 && isequal(t.sleep_s, info.sleepSeconds(1:T)) ...
            && r.S + r.Q + r.W + r.missing == T && r.S_light + r.Q_light + r.W_light + r.missing_light == nnz(t.phase == 1) ...
            && r.S_dark + r.Q_dark + r.W_dark + r.missing_dark == nnz(t.phase == 2) && isequaln(r.r_light, refs.light) && isequaln(r.r_dark, refs.dark) ...
            && all(ismember(t.symbol(~isnan(t.symbol)), [0 1 2])) && r.S == nnz(y == 3);
    end
    ok = ok && string(X.outDir) == string(ld) + "_symbols" && all(ismember([0 1 2], cell2mat(X.symbols)));
    X2 = sqwBatchSymbols(ld); ok = ok && height(X2.summary) == 2;
    if ~ok, why = "a column or count differs from SQWSYMBOLS, or the output folder / rerun is wrong"; end
    [nP, nF] = rep(ok, nP, nF, 2, 'sqwBatchSymbols LD: files = sqwToLabCode(sqwSymbols); counts and references; output beside input', why);
catch e, [nP, nF] = rep(false, nP, nF, 2, 'sqwBatchSymbols LD', e.message); end

%% 3. DD
dd = fullfile(tmp, 'dd'); mkdir(dd); writematrix(syntheticCounts(rs, 2), fullfile(dd, "ddfly.txt"));
try
    X = sqwBatchSymbols(dd, LightDark=false); t = readtable(fullfile(X.outDir, "ddfly_symbols.csv"));
    [~, refs] = sqwSymbols(readCounts(fullfile(dd, "ddfly.txt")), LightDark=false); v = X.summary.Properties.VariableNames;
    ok = all(t.phase == 1) && isequaln(X.summary.r_whole, refs.whole) && ~any(contains(v, ["light" "dark"]));
    [nP, nF] = rep(ok, nP, nF, 3, 'sqwBatchSymbols DD: one phase, r_whole, no light / dark columns', strjoin(string(v), ' '));
catch e, [nP, nF] = rep(false, nP, nF, 3, 'sqwBatchSymbols DD', e.message); end

%% 4. refusals
ids = strings(1, 2);
try, sqwBatchSymbols(fullfile(ld, "fly1.txt")); catch e, ids(1) = e.identifier; end
hd = fullfile(tmp, 'hd'); mkdir(hd); writelines(["count"; string(syntheticCounts(rs, 1))], fullfile(hd, "h.txt"));
try, sqwBatchSymbols(hd); catch e, ids(2) = e.identifier; end
[nP, nF] = rep(ids(1) == "sqw:outDir" && strlength(ids(2)) > 0, nP, nF, 4, 'file list without OutDir and a header line refused', strjoin(ids, ' | '));

%% 5-8. sqwBatchHMM on the LD flies
try
    Rs = sqwBatchHMM(ld + "_symbols", fast{:}, Name="viaSymbols");
    Rc = sqwBatchHMM(ld, fast{:}, Name="viaCounts", OutDir=fullfile(tmp, 'hc'));
    a = Rs.perFly{:, ["ll_K2" "ll_K3" "ll_K4"]}; b = Rc.perFly{:, ["ll_K2" "ll_K3" "ll_K4"]};
    [nP, nF] = rep(isequal(a, b) && isequal(Rs.perFly.file, Rc.perFly.file), nP, nF, 5, ...
        'sqwBatchHMM: symbol files and count files give identical log-likelihoods', mat2str(a - b, 3));
catch e, [nP, nF] = rep(false, nP, nF, 5, 'sqwBatchHMM two routes', e.message); Rs = []; end
if ~isempty(Rs)
    T = Rs.perFly; ok = true;
    for i = 1:height(T)
        ts = readtable(fullfile(ld + "_symbols", T.file(i) + "_symbols.csv")); y = labCodeToSqw(ts.symbol);
        bT = arrayfun(@(K) hmmBic(T.("ll_K" + K)(i), K, 2, numel(y)), [2 3 4]);
        bC = arrayfun(@(K) hmmBic(T.("ll_K" + K)(i), K, 2, [], y), [2 3 4]);
        [~, kT] = min(bT); [~, kC] = min(bC);
        ok = ok && max(abs(bT - T{i, ["BIC_T_K2" "BIC_T_K3" "BIC_T_K4"]})) < 1e-9 && max(abs(bC - T{i, ["BIC_chg_K2" "BIC_chg_K3" "BIC_chg_K4"]})) < 1e-9 ...
            && T.bestK_T(i) == kT + 1 && T.bestK_chg(i) == kC + 1 && T.nchg(i) == nnz(diff(y) ~= 0) ...
            && T.flag_T(i) == any(bT([1 3]) < bT(2) - 10) && T.flag_chg(i) == any(bC([1 3]) < bC(2) - 10) ...
            && abs(T.dBIC_T_K2m3(i) - (bT(1) - bT(2))) < 1e-9 && abs(T.dBIC_T_K4m3(i) - (bT(3) - bT(2))) < 1e-9;
    end
    [nP, nF] = rep(ok, nP, nF, 6, 'BIC (both n), bestK, flag rule and dBIC columns as defined', 'a BIC-derived column disagrees');
    ok = true;
    for i = 1:height(T)
        B = reshape(T{i, "B_H" + [1 1 1 2 2 2 3 3 3] + "_" + ["Q" "W" "S" "Q" "W" "S" "Q" "W" "S"]}, 3, 3)';
        ok = ok && max(abs(sum(B, 2) - 1)) < 1e-9 && B(3, 3) == max(B(:, 3)) && B(2, 2) >= B(1, 2) && T.deg3(i) == isDegenerateH1(B);
        for ph = ["light" "dark"]
            P = reshape(T{i, "P_" + ph + "_H" + repelem(1:3, 3) + "_H" + repmat(1:3, 1, 3)}, 3, 3)';
            ok = ok && max(abs(sum(P, 2) - 1)) < 1e-9;
        end
        ok = ok && isequal(Rs.K3{i}.B, B);
    end
    [nP, nF] = rep(ok, nP, nF, 7, 'labelled K = 3: rows sum to 1, H3 most S, H2 more W, deg3 = isDegenerateH1', 'a labelled row is off');
    ok = true; ord = ["S" "Q" "W"];
    for i = 1:height(T)
        d2 = split(T.dom_K2(i), "+"); d4 = split(T.dom_K4(i), "+");
        [~, o4] = ismember(d4, ord);
        ok = ok && numel(d2) == 2 && numel(d4) == 4 && issorted(o4) && all(o4 > 0) && T.nS_K4(i) + T.nQ_K4(i) + T.nW_K4(i) == 4 ...
            && T.nS_K4(i) == nnz(d4 == "S");
        [~, mx] = max(Rs.fits{i, 3}.B, [], 2); ok = ok && T.nS_K4(i) == nnz(mx == 3);
    end
    [nP, nF] = rep(ok, nP, nF, 8, 'dom_K2 / dom_K4 one symbol per state, sorted; counts add to 4', 'a pattern is malformed');
end

%% 9. resume
try
    d = fullfile(tmp, 'res'); mkdir(d); copyfile(fullfile(ld + "_symbols", "fly1_symbols.csv"), d);
    R1 = sqwBatchHMM(d, fast{:}, Ks=[2 3], Name="r");
    cf = fullfile(d, "r_fits", "fly1.mat"); C = load(cf); C.F{1}.ll = -12345.678; save(cf, '-struct', 'C');
    R2 = sqwBatchHMM(d, fast{:}, Ks=[2 3], Name="r");
    R3 = sqwBatchHMM(d, fast{:}, Ks=[2 3], Name="r", Starts=5);
    ok = R2.perFly.ll_K2 == -12345.678 && abs(R3.perFly.ll_K2 - R1.perFly.ll_K2) < 1 && R3.perFly.ll_K2 ~= -12345.678;   % refitted: near the old optimum, not the planted value
    [nP, nF] = rep(ok, nP, nF, 9, 'Resume: matching cache reused, cache with other Starts ignored', ...
        sprintf('%.3f %.3f %.3f', R1.perFly.ll_K2, R2.perFly.ll_K2, R3.perFly.ll_K2));
catch e, [nP, nF] = rep(false, nP, nF, 9, 'Resume', e.message); end

%% 10. malformed symbol files
cnt = 0; m = (1:1440)';
for c = 1:3
    d = fullfile(tmp, "bad" + c); mkdir(d); s = mod(m, 3); ph = 1 + (m > 720);
    if c == 1, s(5) = 3; elseif c == 2, ph(9) = 3; end
    t = table(m, ph, s, 'VariableNames', {'minute', 'phase', 'symbol'});
    if c == 3, t.phase = []; end
    writetable(t, fullfile(d, "b_symbols.csv"));
    try, sqwBatchHMM(d, fast{:}); catch, cnt = cnt + 1; end
end
[nP, nF] = rep(cnt == 3, nP, nF, 10, 'malformed symbol files refused (code 3, phase 3, no phase column)', sprintf('%d of 3 refused', cnt));

%% 11-12. recovery
for K = [3 2]
    k = 14 - K;
    try
        T = 14400; ph = minutePhases(T);
        if K == 3
            B = [0.85 0.15 0; 0.1 0.9 0; 0 0 1];
            P = cat(3, [0.80 0.08 0.12; 0.02 0.95 0.03; 0.06 0.02 0.92], [0.75 0.05 0.20; 0.02 0.96 0.02; 0.04 0.01 0.95]);
        else
            B = [0.3 0.7 0; 0 0 1]; P = cat(3, [0.95 0.05; 0.04 0.96], [0.97 0.03; 0.02 0.98]);
        end
        p0 = ones(K, 1) / K; [~, y] = hmmSimulatePhased(p0, P, B, ph, RandStream('twister', 'Seed', 40 + K));
        d = fullfile(tmp, "rec" + K); mkdir(d);
        writetable(table((1:T)', ph, sqwToLabCode(y), 'VariableNames', {'minute', 'phase', 'symbol'}), fullfile(d, "sim_symbols.csv"));
        R = sqwBatchHMM(d, Starts=10, Plot=false);
        [~, ~, llTrue] = hmmForward(p0, P, B, y, ph);
        ok = R.perFly.bestK_T == K && R.perFly.("ll_K" + K) >= llTrue - 1e-6;
        [nP, nF] = rep(ok, nP, nF, k, sprintf('recovery: K = %d truth -> bestK_T %d, ll >= ll(true)', K, K), ...
            sprintf('bestK_T %d, ll %.2f vs true %.2f', R.perFly.bestK_T, R.perFly.("ll_K" + K), llTrue));
    catch e, [nP, nF] = rep(false, nP, nF, k, sprintf('recovery K = %d', K), e.message); end
end

fprintf('\n%d passed, %d failed\n', nP, nF);
if nF > 0, error('test_sqw_batch:failed', '%d test(s) failed.', nF); end
end

function x = syntheticCounts(rs, days)
% 1-s counts: alternating sleep (all zero, 6-40 min), quiet wake (p 0.08) and active wake (p 0.6)
n = days * 86400; x = zeros(n, 1); t = 1;
while t <= n
    switch randi(rs, 3)
        case 1, L = 60 * randi(rs, [6 40]); p = 0;
        case 2, L = 60 * randi(rs, [2 10]); p = 0.08;
        otherwise, L = 60 * randi(rs, [2 10]); p = 0.6;
    end
    e = min(n, t + L - 1); x(t:e) = (rand(rs, e - t + 1, 1) < p) .* randi(rs, 4, e - t + 1, 1); t = e + 1;
end
end

function [nP, nF] = rep(ok, nP, nF, k, what, why)
if ok
    nP = nP + 1; fprintf('[PASS] %2d  %s\n', k, what);
else
    nF = nF + 1; fprintf('[FAIL] %2d  %s\n        %s\n', k, what, why);
end
end
