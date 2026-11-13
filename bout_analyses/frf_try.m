cd('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode');
rs = RandStream('mt19937ar', 'Seed', 1);
e = round(-600*log(rs.rand(180,1))); e = e - mean(e);
tic; T = FitResidualFamilies(e); toc
disp(T)
r = GeneralizedHyperbolic_MLE(e, 'NIG', 'NoProfile'); disp(r.Model)
r = GeneralizedHyperbolic_MLE(e, 'HYP', 'NoProfile'); disp(r.Model)
r = GeneralizedHyperbolic_MLE(e, 'GH', 'NoProfile'); disp(r.Model)
