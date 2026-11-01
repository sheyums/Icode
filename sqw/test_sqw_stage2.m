function test_sqw_stage2()
%TEST_SQW_STAGE2  Stage 2 of the sqw_ref port: the HMM bracket, against the golden files and by stress tests.
%
%   Run by name from the sqw folder:   cd <Icode>/sqw; test_sqw_stage2
%
%   Golden (hmm_audit/golden):
%    1. G3  hmmBic parameter counts (9 cases) and BIC value
%    2. G3  hmmViterbiPhased path, exact;  G9 Viterbi with unknown-wake minutes, exact
%    3. G4  labelStates and isDegenerateH1
%    4. G6  hmmFFBS + posteriorSummaries, 4000 paths on G3: mean route counts within 4 SE,
%           c, Xn, Yn, J, EPR within 4 combined Monte Carlo SE (statistical)
%    5. G3  hmmFitPhased reaches a log-likelihood >= that of G3's parameters after 300 EM steps
%   Stress:
%    6. batched hmmFitPhased = iterated hmmEmStep, start by start (with unknown-wake minutes)
%    7. hmmViterbiPhased = brute force over all K^T paths (T = 7, K = 3, two phases, gaps)
%    8. hmmFFBS marginals = exact forward-backward posteriors (one- and two-minute marginals)
%    9. labelStates: invariant to how the fit's states were ordered; likelihood unchanged
%   10. recovery: simulate a known 3-state HMM (8 days), fit, label: P and B within 0.03
%   11. hmmSimulatePhased frequencies within sampling error; posteriorSummaries J = J of mean counts

here = fileparts(mfilename('fullpath'));
assert(isequal(fileparts(which('hmmFitPhased')), here), 'Run from the sqw folder.');
gd = fullfile(here, '..', 'hmm_audit', 'golden');
G = @(f) jsondecode(fileread(fullfile(gd, f)));
nP = 0; nF = 0;
fprintf('\nRunning stage 2 tests for the sqw port...\n\n');
g3 = G('G3_hmm.json');
y3 = pyToSym(g3.y); ph3 = g3.ph + 1; P3 = pyP(g3.P);

%% 1. parameter counts and BIC
try
    ok = true;
    for K = 2:4
        for n = [1 2 4]
            [~, k] = hmmBic(0, K, n, 10);
            ok = ok && k == g3.n_params.(sprintf('K%d_nP%d', K, n));
        end
    end
    b = g3.bic; v = hmmBic(b.ll, b.K, b.nP, b.n);
    ok = ok && abs(v - b.value) < 1e-9;
    [nP, nF] = rep(ok, nP, nF, 1, sprintf('G3 parameter counts (9) and BIC %.6f', v), 'mismatch');
catch err, [nP, nF] = rep(false, nP, nF, 1, '', err.message); end

%% 2. Viterbi
try
    z = hmmViterbiPhased(g3.pi, P3, g3.B, y3, ph3);
    g9 = G('G9_wake_unknown.json');
    z9 = hmmViterbiPhased(g3.pi, P3, g3.B, pyToSym(g9.y), g9.ph + 1);
    [nP, nF] = rep(isequal(z, g3.viterbi + 1) && isequal(z9, g9.viterbi + 1), nP, nF, 2, ...
        'G3 and G9 Viterbi paths exact (2880 minutes each)', ...
        sprintf('G3 mismatches %d, G9 %d', nnz(z ~= g3.viterbi + 1), nnz(z9 ~= g9.viterbi + 1)));
catch err, [nP, nF] = rep(false, nP, nF, 2, '', err.message); end

%% 3. labelling and degeneracy
try
    g = G('G4_matrices.json'); L = g.label;
    [Pl, Bl, pl] = labelStates(pyP(L.P_in), L.B_in, [0.5 0.2 0.3]');
    ok = max(abs(Pl - pyP(L.P_out)), [], 'all') < 1e-15 && isequal(Bl, L.B_out) && ...
        max(abs(pl - L.pi_out)) < 1e-15 && isDegenerateH1(g.degenerate.B_distinct) == g.degenerate.distinct && ...
        isDegenerateH1(g.degenerate.B_copy) == g.degenerate.copy;
    [nP, nF] = rep(ok, nP, nF, 3, 'G4 labelling (P, B, pi) and degeneracy (distinct / copy)', 'mismatch');
catch err, [nP, nF] = rep(false, nP, nF, 3, '', err.message); end

%% 4. FFBS against G6 (statistical)
try
    g6 = G('G6_ffbs.json');
    Z = hmmFFBS(g3.pi, P3, g3.B, y3, ph3, NPaths=4000, RandomSeed=6);
    hrs = accumarray(ph3, 1, [2 1]) / 60;
    S = posteriorSummaries(Z, ph3, 2, hrs);
    rcK = zeros(4000, 2, 4);
    for k = 1:4000, rcK(k, :, :) = routeCounts(Z(k, :)', ph3, 2); end
    se = squeeze(std(rcK, 0, 1)) / sqrt(4000);
    zr = (S.routeCounts - g6.mean_route_counts) ./ hypot(se, g6.se_route_counts);
    % per-path statistics: Monte Carlo SE from the port's own paths, doubled for the Python's
    Js = zeros(4000, 2); Es = Js;
    for k = 1:4000
        N = pairCounts(Z(k, :)', ph3, 2); Js(k, :) = jCycleCount(N, hrs)';
        for p = 1:2, e = epr_count_based(N(:, :, p), hrs(p) * 3600, 'bound'); Es(k, p) = e; end
    end
    zJ = (S.J - g6.J) ./ (sqrt(2) * std(Js)' / sqrt(4000));
    zE = (S.EPR_bound - g6.EPR_bound) ./ (sqrt(2) * std(Es)' / sqrt(4000));
    dXYc = max(abs([S.Xn - g6.Xn; S.Yn - g6.Yn; S.c - g6.c]));
    ok = all(abs(zr(:)) < 4) && all(abs(zJ) < 4) && all(abs(zE) < 4) && dXYc < 0.01;
    [nP, nF] = rep(ok, nP, nF, 4, sprintf(['G6 FFBS 4000 paths: route counts max |z| %.2f, J |z| %.2f, ' ...
        'EPR |z| %.2f, Xn/Yn/c max diff %.4f'], max(abs(zr(:))), max(abs(zJ)), max(abs(zE)), dXYc), ...
        sprintf('route z %s, J z %s, EPR z %s, XYc %.4f', mat2str(zr, 2), mat2str(zJ', 2), mat2str(zE', 2), dXYc));
catch err, [nP, nF] = rep(false, nP, nF, 4, '', err.message); end

%% 5. the fit reaches at least the optimum near G3's parameters
try
    [p1, Pk, Bk] = deal(g3.pi, P3, g3.B);
    for it = 1:300, [p1, Pk, Bk, llRef] = hmmEmStep(p1, Pk, Bk, y3, ph3); end
    fit = hmmFitPhased(y3, ph3, 2, 3, Starts=20, RandomSeed=20);
    [nP, nF] = rep(fit.ll >= llRef - 1e-6, nP, nF, 5, sprintf(['hmmFitPhased on G3: ll %.6f >= %.6f ' ...
        '(G3 parameters after 300 EM steps); %d iterations'], fit.ll, llRef, fit.iterations), ...
        sprintf('fit %.6f vs reference %.6f', fit.ll, llRef));
catch err, [nP, nF] = rep(false, nP, nF, 5, '', err.message); end

%% 6. batched EM = iterated single steps
try
    g9 = G('G9_wake_unknown.json'); y9 = pyToSym(g9.y); ph9 = g9.ph + 1;
    rs = RandStream('twister', 'Seed', 61); dmax = 0;
    for r = 1:3
        st.pi = randDir(rs, 3, 1); st.P = permute(reshape(randDir(rs, 3, 6), 3, 3, 2), [2 1 3]);
        st.B = randDir(rs, 3, 3)';
        f = hmmFitPhased(y9, ph9, 2, 3, Start=st, MaxIter=40, Tol=1e-30);
        [a, b, c] = deal(st.pi, st.P, st.B);
        for it = 1:40, [a, b, c] = hmmEmStep(a, b, c, y9, ph9); end
        dmax = max([dmax, max(abs(f.P - b), [], 'all'), max(abs(f.B - c), [], 'all'), max(abs(f.pi - a))]);
    end
    [nP, nF] = rep(dmax < 1e-10, nP, nF, 6, sprintf('batched fit = 40 iterated EM steps, 3 random starts, unknown wake (max diff %.1e)', dmax), ...
        sprintf('max diff %.3g', dmax));
catch err, [nP, nF] = rep(false, nP, nF, 6, '', err.message); end

%% 7. Viterbi vs brute force
try
    rs = RandStream('twister', 'Seed', 7); ok = true;
    for r = 1:50
        T = 7; pi0 = randDir(rs, 3, 1); P = permute(reshape(randDir(rs, 3, 6), 3, 3, 2), [2 1 3]);
        B = randDir(rs, 3, 3)'; y = randi(rs, [0 4], T, 1); ph = [1 1 1 2 2 2 1]';
        E = hmmEmission(B, y); best = -Inf;
        for code = 0:3 ^ T - 1
            z = mod(floor(code ./ 3 .^ (T - 1:-1:0)'), 3) + 1;
            lp = log(pi0(z(1))) + log(E(1, z(1)));
            for t = 2:T, lp = lp + log(P(z(t - 1), z(t), ph(t))) + log(E(t, z(t))); end
            if lp > best + 1e-12, best = lp; zb = z; end
        end
        ok = ok && isequal(hmmViterbiPhased(pi0, P, B, y, ph), zb);
    end
    [nP, nF] = rep(ok, nP, nF, 7, 'Viterbi = brute force over 3^7 paths (50 random models, gaps, unknown wake)', 'mismatch');
catch err, [nP, nF] = rep(false, nP, nF, 7, '', err.message); end

%% 8. FFBS marginals vs exact posteriors
try
    rs = RandStream('twister', 'Seed', 8);
    T = 60; pi0 = [0.3 0.3 0.4]'; P = cat(3, [0.8 0.1 0.1; 0.2 0.7 0.1; 0.1 0.1 0.8], [0.6 0.3 0.1; 0.1 0.8 0.1; 0.2 0.2 0.6]);
    B = [0.7 0.25 0.05; 0.2 0.75 0.05; 0.05 0.05 0.9]; ph = 1 + (mod(0:T - 1, 20) >= 10)';
    y = randi(rs, 3, T, 1); y(20:24) = 0; y(40) = 4;
    [A, c] = hmmForward(pi0, P, B, y, ph); E = hmmEmission(B, y);
    Bk = ones(T, 3);
    for t = T - 1:-1:1, Bk(t, :) = ((P(:, :, ph(t + 1)) * (E(t + 1, :) .* Bk(t + 1, :))') / c(t + 1))'; end
    Gm = A .* Bk; Gm = Gm ./ sum(Gm, 2);
    t2 = 30; X2 = (A(t2, :)' * (E(t2 + 1, :) .* Bk(t2 + 1, :))) .* P(:, :, ph(t2 + 1)); X2 = X2 / sum(X2(:));
    n = 20000; Z = hmmFFBS(pi0, P, B, y, ph, NPaths=n, RandomSeed=8);
    Gh = zeros(T, 3); for k = 1:3, Gh(:, k) = mean(Z == k, 1)'; end
    X2h = accumarray([Z(:, t2), Z(:, t2 + 1)], 1, [3 3]) / n;
    zG = max(abs(Gh - Gm) ./ sqrt(max(Gm .* (1 - Gm), 1e-12) / n), [], 'all');
    zX = max(abs(X2h - X2) ./ sqrt(max(X2 .* (1 - X2), 1e-12) / n), [], 'all');
    [nP, nF] = rep(zG < 5 && zX < 5, nP, nF, 8, sprintf(['FFBS marginals = exact posteriors (20000 paths): ' ...
        'max |z| one-minute %.2f, two-minute %.2f'], zG, zX), sprintf('zG %.2f zX %.2f', zG, zX));
catch err, [nP, nF] = rep(false, nP, nF, 8, '', err.message); end

%% 9. labelling invariance
try
    rs = RandStream('twister', 'Seed', 9); ok = true;
    B0 = [0.75 0.2 0.05; 0.1 0.85 0.05; 0.03 0.02 0.95];
    P0 = cat(3, [0.8 0.15 0.05; 0.06 0.9 0.04; 0.02 0.03 0.95], [0.7 0.2 0.1; 0.1 0.85 0.05; 0.03 0.05 0.92]);
    pi0 = [0.2 0.3 0.5]';
    ll0 = llOf(pi0, P0, B0, y3, ph3);
    for r = 1:20
        o = randperm(rs, 3);
        [Pl, Bl, pl] = labelStates(P0(o, o, :), B0(o, :), pi0(o));
        ok = ok && isequal(Pl, P0) && isequal(Bl, B0) && isequal(pl, pi0) && ...
            abs(llOf(pi0(o), P0(o, o, :), B0(o, :), y3, ph3) - ll0) < 1e-9;
    end
    [nP, nF] = rep(ok, nP, nF, 9, 'labelStates recovers (H1, H2, H3) from 20 random orderings; likelihood invariant', 'mismatch');
catch err, [nP, nF] = rep(false, nP, nF, 9, '', err.message); end

%% 10. recovery of a known HMM
try
    rs = RandStream('twister', 'Seed', 10);
    ph = minutePhases(8 * 1440);
    [~, y] = hmmSimulatePhased(pi0, P0, B0, ph, rs);
    f = hmmFitPhased(y, ph, 2, 3, Starts=10, RandomSeed=1);
    [Pf, Bf] = labelStates(f.P, f.B, f.pi);
    dP = max(abs(Pf - P0), [], 'all'); dB = max(abs(Bf - B0), [], 'all');
    llTrue = llOf(pi0, P0, B0, y, ph);
    [nP, nF] = rep(dP < 0.03 && dB < 0.03 && f.ll >= llTrue - 1e-6, nP, nF, 10, sprintf(['recovery on 8 simulated days: ' ...
        'max |dP| %.3f, |dB| %.3f; fitted ll %.2f >= true-parameter ll %.2f'], dP, dB, f.ll, llTrue), ...
        sprintf('dP %.3f dB %.3f ll %.2f vs %.2f', dP, dB, f.ll, llTrue));
catch err, [nP, nF] = rep(false, nP, nF, 10, '', err.message); end

%% 11. simulation frequencies; J linear in counts
try
    rs = RandStream('twister', 'Seed', 11);
    ph = ones(200000, 1); P1 = P0(:, :, 1);
    [z, y] = hmmSimulatePhased(pi0, P1, B0, ph, rs);
    Nt = accumarray([z(1:end - 1), z(2:end)], 1, [3 3]); Pe = Nt ./ sum(Nt, 2);
    Ne = accumarray([z, y], 1, [3 3]); Be = Ne ./ sum(Ne, 2);
    zP = max(abs(Pe - P1) ./ sqrt(P1 .* (1 - P1) ./ sum(Nt, 2)), [], 'all');
    zB = max(abs(Be - B0) ./ sqrt(B0 .* (1 - B0) ./ sum(Ne, 2)), [], 'all');
    Z = hmmFFBS(g3.pi, P3, g3.B, y3, ph3, NPaths=50, RandomSeed=1);
    hrs = accumarray(ph3, 1, [2 1]) / 60;
    S = posteriorSummaries(Z, ph3, 2, hrs);
    okJ = max(abs(S.J - jCycleCount(S.pairCounts, hrs))) < 1e-12;
    [nP, nF] = rep(zP < 5 && zB < 5 && okJ, nP, nF, 11, sprintf(['simulation: transition / emission frequencies ' ...
        'max |z| %.2f / %.2f; posterior J = J of mean counts'], zP, zB), sprintf('zP %.2f zB %.2f J %d', zP, zB, okJ));
catch err, [nP, nF] = rep(false, nP, nF, 11, '', err.message); end

fprintf('\n%d passed, %d failed\n', nP, nF);
if nF > 0, error('test_sqw_stage2:failed', '%d test(s) failed', nF); end
end

% ----------------------------------------------------------------------- helpers
function [nP, nF] = rep(ok, nP, nF, k, what, why)
if ok
    nP = nP + 1; fprintf('[PASS] %2d  %s\n', k, what);
else
    nF = nF + 1; fprintf('[FAIL] %2d  %s\n        %s\n', k, what, why);
end
end
function y = pyToSym(v), y = v(:) + 1; y(v(:) == -1) = 0; y(v(:) == -2) = 4; end
function P = pyP(v), P = permute(v, [2 3 1]); end
function X = randDir(rs, k, n), g = -log(rand(rs, k, n)); X = g ./ sum(g, 1); end
function ll = llOf(pi0, P, B, y, ph), [~, ~, ll] = hmmForward(pi0, P, B, y, ph); end
