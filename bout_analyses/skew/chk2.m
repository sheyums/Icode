addpath('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode');
rs = RandStream('mt19937ar','Seed',1);
P = struct('mu',0,'alpha',2,'beta',0.8,'delta',1); P.gamma = sqrt(P.alpha^2-P.beta^2);
m=P.delta/P.gamma; lam=P.delta^2; y=randn(rs,800,1).^2;
v=m+(m^2*y)/(2*lam)-(m/(2*lam))*sqrt(4*m*lam*y+m^2*y.^2); u=rand(rs,800,1); fl=u>m./(m+v); v(fl)=m^2./v(fl);
x=P.mu+P.beta*v+sqrt(v).*randn(rs,800,1);
for c = [1 1e-6]
 r = GeneralizedHyperbolic_MLE(c*x-7*(c~=1),'NIG','NoProfile');
 fprintf('c=%g SEValid %d SE %s skew %.6g skSE %.4g deltaAtBound %d\n', c, r.SEValid, mat2str(r.SE',4), r.Skewness, r.SkewnessSE, r.delta_at_bound);
end
