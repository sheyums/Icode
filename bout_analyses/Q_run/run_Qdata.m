[X, ds] = qload(); R = qrun(X, 1000, 6);
[kap, ci, kb] = qkappa(R.fast, ds, 2000, 7);
dn = ["per0 DD from 24 h", "CS LD", "CS DD from 24 h"];
fprintf('#### pooled %d flies: fast F1 %.4f F2 %.4f F3 %.4f | sleep-sleep %.4f wake-wake %.4f\n', numel(X), R.F(1:5));
fprintf('  raw r1 %.4f r2 %.4f r3 %.4f; W2-null means %.4f %.4f %.4f\n', mean(R.obs(:, 1:3), 'omitnan'), mean(R.nullmean(:, 1:3), 'omitnan'));
fprintf('  T1 (F2 > 0): mean r2 %.4f vs W2 null %.4f, p %.4f\n', mean(R.obs(isfinite(R.obs(:, 2)), 2)), mean(R.g2null), R.pT1);
fprintf('  T2 kappa = F2/F1 = %.3f, 95%% CI %.3f..%.3f -> ', kap, ci);
if ci(2) < 0.5, fprintf('PAIR-LOCAL\n'); elseif ci(1) > 0.5, fprintf('MULTI-BOUT\n'); else, fprintf('UNRESOLVED\n'); end
k3 = mean(R.fast(:, 3), 'omitnan') / mean(R.fast(:, 1), 'omitnan'); fprintf('  kappa3 = F3/F1 = %.3f\n', k3);
for k = 1:3
    m = ds == k; kk = mean(R.fast(m, 2), 'omitnan') / mean(R.fast(m, 1), 'omitnan');
    fprintf('  %-18s F1 %.4f F2 %.4f F3 %.4f SS %.4f WW %.4f kappa %.3f\n', dn(k), mean(R.fast(m, 1:5), 'omitnan'), kk);
end
save('Q_results.mat', 'R', 'ds', 'kap', 'ci', 'kb');
