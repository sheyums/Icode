function O = o1c_run(csvFile, outDir)
%O1C_RUN  PRESPEC O1-C (beam_hmm/PRESPEC.md, 2026-11-19): within fly x light phase x bout type, one exponential vs a
%two-component hyperexponential, parametric bootstrap LRT (repo-root HyperexponentialLRT, B = 499).
%
%   Run from the Icode REPO ROOT in MATLAB R2025a+ (not Octave):
%     addpath('beam_hmm');
%     O = o1c_run('beam_hmm/results/o1a/bouts_weeks2-3.csv', 'beam_hmm/results/o1c');
%   Writes outDir/o1c_report.txt (full Command Window output), outDir/o1c_units.csv (one row per unit),
%   outDir/o1c_results.mat.

assert(isequal(fileparts(which('HyperexponentialLRT')), pwd), ...
       'Run from the Icode repo root, so the repo copies of the bout functions are used (CLAUDE.md).');
if ~exist(outDir, 'dir'), mkdir(outDir); end
logFile = fullfile(outDir, 'o1c_report.txt'); if exist(logFile, 'file'), delete(logFile); end
diary(logFile); cleanup = onCleanup(@() diary('off'));
T = readtable(csvFile, 'TextType', 'string');
fprintf('O1-C  %s  (%s)\n', csvFile, char(datetime('now')));
kinds = ["sleep", "wake"]; xmins = [300 2]; phases = ["light", "dark"]; flies = unique(T.ch)';
rows = {}; u = 0; untested = 0;
for ik = 1:2
    for c = flies
        for ip = 1:2
            u = u + 1;
            m = T.kind == kinds(ik) & T.ch == c & ((T.start_zt < 12) == (ip == 1));
            x = double(T.length_s(m)); n = numel(x);
            r = struct('kind', kinds(ik), 'ch', c, 'phase', phases(ip), 'n', n, 'LR', NaN, 'p', NaN, 'AtomAtZero', NaN, ...
                       'tau1', NaN, 'tau2', NaN, 'q1', NaN, 'q2', NaN, 'tau1_K1', NaN, 'p23', NaN, 'BValid', NaN);
            if n < 50
                untested = untested + 1; rows{end + 1} = r; %#ok<AGROW>
                fprintf('%-5s ch%d %-5s n %4d  untested (< 50 bouts)\n', kinds(ik), c, phases(ip), n); continue
            end
            L = HyperexponentialLRT(x, xmins(ik), 1, 2, SamplingInterval=1, B=499, RandomSeed=20261120 + u, Verbose=false);
            r.LR = L.LR; r.p = L.pValue; r.AtomAtZero = L.AtomAtZero; r.BValid = L.BValid;
            r.tau1_K1 = L.Fit0.Tau(1);
            [tau, o] = sort(L.Fit1.Tau(:)'); q = L.Fit1.WeightsObserved(:)'; q = q(o);
            r.tau1 = tau(1); r.tau2 = tau(end); r.q1 = q(1); r.q2 = q(end);
            if n >= 150 && r.p < 0.05                                  % reported only
                L3 = HyperexponentialLRT(x, xmins(ik), 2, 3, SamplingInterval=1, B=499, RandomSeed=20262120 + u, Verbose=false);
                r.p23 = L3.pValue;
            end
            fprintf(['%-5s ch%d %-5s n %4d  LR %8.2f  p %.3f  (valid B %d, atom %.2f) | K1 tau %7.1f s | K2 tau %7.1f / ' ...
                     '%7.1f s, q %.2f / %.2f | p(2 vs 3) %s\n'], kinds(ik), c, phases(ip), n, r.LR, r.p, r.BValid, ...
                    r.AtomAtZero, r.tau1_K1, r.tau1, r.tau2, r.q1, r.q2, num2str(r.p23, 3));
            rows{end + 1} = r; %#ok<AGROW>
        end
    end
end
U = struct2table([rows{:}]);
writetable(U, fullfile(outDir, 'o1c_units.csv'));
fprintf('\nuntested units (< 50 bouts): %d\n', untested);
for ik = 1:2
    s = U(U.kind == kinds(ik) & ~isnan(U.p), :); mT = height(s); rr = nnz(s.p < 0.05);
    if rr == 0, pb = 1; else, pb = betainc(0.05, rr, mT - rr + 1); end   % P(X >= r), X ~ Bin(m, 0.05); no toolbox
    if pb < 0.05, v = 'SUPPORTED'; else, v = 'NOT SUPPORTED'; end
    fprintf('>> %s: %d of %d fly-phases reject one time scale (p < 0.05); binomial p (vs 0.05) = %.3g -> within-fly-phase multiple time scales %s\n', ...
            upper(kinds(ik)), rr, mT, pb, v);
    for ip = 1:2
        sp = s(s.phase == phases(ip), :);
        fprintf('   %s: %d of %d reject; K2 tau1 median %.1f s (range %.1f-%.1f), tau2 median %.1f s (range %.1f-%.1f), q1 median %.2f\n', ...
                phases(ip), nnz(sp.p < 0.05), height(sp), median(sp.tau1), min(sp.tau1), max(sp.tau1), ...
                median(sp.tau2), min(sp.tau2), max(sp.tau2), median(sp.q1));
    end
    O.(char(kinds(ik))) = struct('m', mT, 'r', rr, 'binomial_p', pb, 'verdict', v);
end
O.units = U;
save(fullfile(outDir, 'o1c_results.mat'), 'O');
fprintf('\nsaved %s\n', fullfile(outDir, 'o1c_results.mat'));
end
