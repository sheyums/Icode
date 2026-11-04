function S = sleepSensitivityTests(Tobs, Thid, options)
%SLEEPSENSITIVITYTESTS  Post hoc 22h sensitivity: hidden J against OBSERVED sleep, removing the shared HMM fit from x (ledger 2026-11-03).
%
%   S = sleepSensitivityTests(Tobs, Thid)
%   S = sleepSensitivityTests(Tobs, Thid, NPerm=10000, Seed=223)
%
%   Tobs   SLEEPVSCYCLEROWS(D) (observed sleep)
%   Thid   SLEEPVSCYCLEROWS(D, Level="hidden", ...) (model sleep 60 pi_S, J_hmm, cycle)
%   Rows are matched by genotype, fly and phase. Per PRESPEC 22 primary block:
%     rho_sleep   Spearman of model sleep against observed sleep (all matched rows)
%     rho_J, p_J  Spearman of J_hmm against OBSERVED sleep (rows with a hidden cycle), permutation p
%                 (SPEARMANPERMTEST, values rounded to 1e-10 before ranking)
%   Descriptive, post hoc; no Holm. Expected (ledger): model vs observed sleep 0.99-1.00; J against observed
%   sleep sss L -0.49 (p 0.003), fumin D +0.53, fumin L +0.40, sss D +0.14 -- the model-sleep values.
%   S  table: genotype, phase, n, rho_sleep, n_J, rho_J, p_J
%
%   Example (from the Icode\sqw folder):
%       To = sleepVsCycleRows(D); Th = sleepVsCycleRows(D, Level="hidden", DataRoot="D:\flies", JsonlFile=jsonl);
%       S = sleepSensitivityTests(To, Th)
%
%   See also SLEEPVSCYCLEROWS, SLEEPVSCYCLETESTS, SPEARMANPERMTEST.

arguments
    Tobs table
    Thid table
    options.NPerm (1,1) double = 10000
    options.Seed (1,1) double = 223
end
rs = RandStream('twister', 'Seed', options.Seed);
o = Tobs(:, {'genotype', 'fly', 'phase', 'sleep'}); o.Properties.VariableNames{'sleep'} = 'sleepObs';
M = innerjoin(Thid, o, 'Keys', {'genotype', 'fly', 'phase'});
prim = {"CS", 1; "CS", 2; "per0DD", 1; "fumin", 1; "fumin", 2; "sss", 1; "sss", 2}; S = table();
for k = 1:size(prim, 1)
    U = M(M.genotype == prim{k, 1} & M.phase == prim{k, 2}, :); C = U(U.cycle, :);
    [rj, pj] = spearmanPermTest(C.sleepObs, C.J, options.NPerm, rs);
    S = [S; table(prim{k, 1}, prim{k, 2}, height(U), spearmanRho(round(U.sleep, 10), round(U.sleepObs, 10)), height(C), rj, pj, ...
        'VariableNames', {'genotype', 'phase', 'n', 'rho_sleep', 'n_J', 'rho_J', 'p_J'})]; %#ok<AGROW>
end
end
