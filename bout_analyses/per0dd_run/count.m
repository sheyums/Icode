out = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\per0dd_run';
for t = ["sleep","wake"]
  T = readtable(fullfile(out, "per0dd_" + t + "_perfly.csv"));
  fprintf('%s:\n', t);
  for k = [1 2 4 8]
    r = T(T.Lag == k, :);
    fprintf('  lag %d: p<0.05 IAAFT  Up %2d (up>0.5: %d, <0.5: %d)  Bowley %2d | Rho shuffle %2d (rho>0: %d)  | Up range %.3f-%.3f, median NPairs %g\n', k, ...
      sum(r.UpFraction_P_IAAFT < 0.05), sum(r.UpFraction_P_IAAFT < 0.05 & r.UpFraction > 0.5), sum(r.UpFraction_P_IAAFT < 0.05 & r.UpFraction < 0.5), ...
      sum(r.Bowley_P_IAAFT < 0.05), sum(r.Rho_P_Shuffle < 0.05), sum(r.Rho_P_Shuffle < 0.05 & r.Rho > 0), min(r.UpFraction), max(r.UpFraction), median(r.NPairs));
  end
end
load(fullfile(out, 'per0dd_results.mat'), 'res');
for t = ["sleep","wake"]
  B = res.(t + "_bouts"); d = diff(B, 1, 2); d = d(~isnan(d));
  fprintf('%s lag-1 differences: %d, zero %d (%.1f%%)\n', t, numel(d), sum(d == 0), 100*mean(d == 0));
end
