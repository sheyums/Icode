% set SC (1-3) and NG before calling
[X, ds] = qload(); nm = ["a independent", "b shared fast state rho 0.8", "c one-step chain rho 0.06"];
rej = 0; cls = zeros(1, 3); K = NaN(NG, 1); F12 = NaN(NG, 2);
for g = 1:NG
    rs = RandStream('mrg32k3a', 'Seed', 60000 + 1000 * SC + g);
    Y = cellfun(@(x) qsim(x, SC, rs), X, 'UniformOutput', false);
    R = qrun(Y, 199, 70000 + 1000 * SC + g);
    [K(g), ci] = qkappa(R.fast, ds, 500, 80000 + 1000 * SC + g);
    rej = rej + (R.pT1 < 0.05); F12(g, :) = R.F(1:2);
    if ci(2) < 0.5, cls(1) = cls(1) + 1; elseif ci(1) > 0.5, cls(2) = cls(2) + 1; else, cls(3) = cls(3) + 1; end
end
fprintf('CAL %-30s T1 rejection %.2f | kappa: pair-local %.2f, multi-bout %.2f, unresolved %.2f | median kappa %.3f | mean F1 %.4f F2 %.4f\n', ...
    nm(SC), rej / NG, cls / NG, median(K), mean(F12));
save(sprintf('Qcal%d.mat', SC), 'rej', 'cls', 'K', 'F12', 'NG');
