cd('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts');
w = which('GeneralizedHyperbolic_MLE'); disp(w)
assert(strcmp(fileparts(w), pwd))
rng(3); x = 0.8*gamrnd(1,1,600,1) + sqrt(gamrnd(1,1,600,1)).*randn(600,1);
r = GeneralizedHyperbolic_MLE(x, 'NIG', 'NoProfile');
fprintf('Skew %.4f SE %.4f BowleyModel %.4f BowleyData %.4f SEValid %d\n', r.Skewness, r.SkewnessSE, r.SkewnessBowleyModel, r.SkewnessBowleyData, r.SEValid);
