function [B, F] = completedLoopRates(D)
%COMPLETEDLOOPRATES  Completed cycles per hour, hours per net cycle and loop durations per genotype x phase (cycle_times.py).
%
%   [B, F] = completedLoopRates(D)
%
%   D   BUILDPERFLYDATA output (six genotypes; per fly the symbols y and phases ph, per fly x phase
%       the Q/W reference flag refOK and the chain's J_chain)
%
%   Per fly x phase WITH a Q/W reference: completed cycles (COMPLETEDLOOPS) assigned to that phase,
%   divided by the scored hours in it (minutes with a symbol, / 60). Descriptive; no test.
%   B   table, one row per genotype x phase (cycle_times.json): genotype, phase, n (flies),
%       fwd_per_h, back_per_h, net_per_h (means over flies), J_chain_mean (mean chain J of the same
%       flies, no-cycle = 0; net_per_h should equal it), h_per_net_cycle (1 / net_per_h if > 0, else
%       NaN), dur_fwd_median_min, dur_back_median_min, dur_fwd_mean_min, dur_back_mean_min (pooled
%       over the flies' loops), n_loops_fwd, n_loops_back
%   F   table, one row per fly x phase used: genotype, fly, phase, hours, fwd, back (counts),
%       fwd_per_h, back_per_h, net_per_h, J_chain
%
%   Example (from the Icode\sqw folder):
%       D = buildPerFlyData("D:\flies", fullfile("..", "hmm_audit", "sqw20b_fit2_theta0.5.jsonl"), Cache="perfly_cache.mat");
%       B = completedLoopRates(D); disp(B(:, ["genotype" "phase" "net_per_h" "J_chain_mean"]))
%
%   See also COMPLETEDLOOPS, CYCLEWALK, EXPORTCUMULATIVEWALKFIGURE.

arguments
    D (1,1) struct
end
B = table(); F = table();
for k = 1:numel(D.G)
    G = D.G(k); T = G.perFly; g = D.order(k);
    Lf = arrayfun(@(f) completedLoops(f.y, f.ph), G.flies, 'UniformOutput', false);
    for p = unique(T.phase)'
        R = T(T.phase == p & T.refOK, :); df = []; db = []; rows = cell(height(R), 1);
        for r = 1:height(R)
            fl = G.flies(R.fly(r)); L = Lf{R.fly(r)}; L = L(L.phase == p, :);
            h = nnz(fl.ph == p & fl.y >= 1 & fl.y <= 3) / 60;
            a = nnz(L.dir == 1); b = nnz(L.dir == -1); J = R.J_chain(r); if ~isfinite(J), J = 0; end
            rows{r} = table(g, R.fly(r), p, h, a, b, a / h, b / h, (a - b) / h, J, 'VariableNames', ...
                {'genotype', 'fly', 'phase', 'hours', 'fwd', 'back', 'fwd_per_h', 'back_per_h', 'net_per_h', 'J_chain'});
            df = [df; L.dur(L.dir == 1)]; db = [db; L.dur(L.dir == -1)]; %#ok<AGROW>
        end
        Fp = vertcat(rows{:}); F = [F; Fp]; %#ok<AGROW>
        net = mean(Fp.net_per_h); hpc = NaN; if net > 0, hpc = 1 / net; end
        B = [B; table(g, p, height(Fp), mean(Fp.fwd_per_h), mean(Fp.back_per_h), net, mean(Fp.J_chain), hpc, ...
            median(df), median(db), mean(df), mean(db), numel(df), numel(db), 'VariableNames', {'genotype', 'phase', 'n', ...
            'fwd_per_h', 'back_per_h', 'net_per_h', 'J_chain_mean', 'h_per_net_cycle', 'dur_fwd_median_min', ...
            'dur_back_median_min', 'dur_fwd_mean_min', 'dur_back_mean_min', 'n_loops_fwd', 'n_loops_back'})]; %#ok<AGROW>
    end
end
end
