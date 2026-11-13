cd('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode');
rs = RandStream('mt19937ar', 'Seed', 1);
for n = [90 180]
  e = round(-600 * log(rs.rand(n, 1))); e = e - mean(e);
  tic; r = GeneralizedHyperbolic_MLE(e, 'AUTO', 'NoProfile'); t = toc;
  fprintf('n %d: AUTO NoProfile %.1f s, chose %s, logL %.2f, skew %.2f\n', n, t, string(r.ModelSelection), r.LogLike, r.Skewness);
end
