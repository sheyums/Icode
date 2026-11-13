out = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\ar1_run';
f = dir(fullfile(out, '*_perfly.csv'));
for i = 1:numel(f)
  T = readtable(fullfile(out, f(i).name)); T = T(isfinite(T.B), :);
  s = T.Hetero_P < 0.05; b = T.B_P < 0.05;
  fprintf('%-26s n %2d | B sig %2d (pos %2d) | Hetero sig %d (pos %d) | ShapeKS sig %2d | BelowHalfB med %.3f (null %.3f) | NPairs med %g\n', ...
    erase(f(i).name, '_perfly.csv'), height(T), sum(b), sum(b & T.B > 0), sum(s), sum(s & T.Hetero > 0), ...
    sum(T.ShapeKS_P < 0.05), median(T.BelowHalfB, 'omitnan'), median(T.BelowHalfB_NullMean, 'omitnan'), median(T.NPairs));
end
load(fullfile(out, 'ar1_results.mat'), 'allRes');
G = allRes.per0dd_sleep_all.R.Group; fprintf('per0 sleep group B mean %.3f p %.3f; ', G.B_Mean, G.B_Mean_P);
G = allRes.per0dd_wake_all.R.Group; fprintf('per0 wake group B mean %.3f p %.3f\n', G.B_Mean, G.B_Mean_P);
for t = string(fieldnames(allRes)).'
  G = allRes.(t).R.Group; fprintf('%-20s B NSig %d CountP %.4f | ShapeKS NSig %d CountP %.4f\n', t, G.B_NSig, G.B_CountP, G.ShapeKS_NSig, G.ShapeKS_CountP);
end
