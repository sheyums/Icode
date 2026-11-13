% Independent checks of the new skewness outputs
addpath('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode');
rs = RandStream('mt19937ar','Seed',1);
% NIG sample
P = struct('mu',0,'alpha',2,'beta',0.8,'delta',1); P.gamma = sqrt(P.alpha^2-P.beta^2);
x = nigs(rs,P,800);
for flag = {'NIG','HYP','GH'}
  r = GeneralizedHyperbolic_MLE(5*x+3, flag{1}, 'NoProfile');
  p = num2cell(r.Params); [mu,lam,al,be,de] = p{:}; ga = sqrt(al^2-be^2);
  % (1) symbolic exact moment skewness
  digits(50); w = sym(de)*sym(ga); L = sym(lam);
  m = @(k) (sym(de)/sym(ga))^k*besselk(L+k,w)/besselk(L,w);
  k1=m(1); k2=m(2)-k1^2; k3=m(3)-3*k1*m(2)+2*k1^3; b=sym(be);
  gs = double(vpa((3*b*k2+b^3*k3)/(k1+b^2*k2)^(3/2)));
  % (2) NIG closed form
  gn = NaN; if lam==-0.5, gn = 3*be/(al*sqrt(de*ga)); end
  % (3) model Bowley by brute trapz on a fine grid
  xs = linspace(mu-60*de-60/ (al-abs(be)), mu+60*de+60/(al-abs(be)), 2e6);
  f = exp(ghlog(xs,mu,lam,al,be,de)); F = cumtrapz(xs,f); F=F/F(end);
  [Fu,iu]=unique(F); qq = interp1(Fu,xs(iu),[.25 .5 .75]); bm = (qq(3)+qq(1)-2*qq(2))/(qq(3)-qq(1));
  % (4) data Bowley vs quantile()
  y = 5*x+3; qd = quantile(y,[.25 .5 .75]); bd = (qd(3)+qd(1)-2*qd(2))/(qd(3)-qd(1));
  fprintf('%s: skew %.10g sym %.10g nig %.10g | SE %.4g CI [%.4g %.4g] | BowleyModel %.8g brute %.8g | BowleyData %.10g quantile %.10g\n', ...
    flag{1}, r.Skewness, gs, gn, r.SkewnessSE, r.SkewnessCIWald, r.SkewnessBowleyModel, bm, r.SkewnessBowleyData, bd);
end
% Equivariance: scale/shift/negation
r1 = GeneralizedHyperbolic_MLE(x,'NIG','NoProfile'); r2 = GeneralizedHyperbolic_MLE(1e-6*x-7,'NIG','NoProfile'); r3 = GeneralizedHyperbolic_MLE(-x,'NIG','NoProfile');
fprintf('equivariance: %.3g %.3g %.3g | %.3g %.3g %.3g | neg: %.3g %.3g %.3g\n', r1.Skewness-r2.Skewness, r1.SkewnessSE-r2.SkewnessSE, r1.SkewnessBowleyModel-r2.SkewnessBowleyModel, ...
  r1.SkewnessBowleyData-r2.SkewnessBowleyData, 0,0, r1.Skewness+r3.Skewness, r1.SkewnessBowleyModel+r3.SkewnessBowleyModel, r1.SkewnessBowleyData+r3.SkewnessBowleyData);
% Large w (near Normal) symbolic check of closed form, via a fit on normal data
z = randn(rs,2000,1)+0.05*randn(rs,2000,1).^2;
r = GeneralizedHyperbolic_MLE(z,'NIG','NoProfile'); p=num2cell(r.Params); [mu,lam,al,be,de]=p{:}; ga=sqrt(al^2-be^2);
w = sym(de)*sym(ga); m = @(k) (sym(de)/sym(ga))^k*besselk(sym(lam)+k,w)/besselk(sym(lam),w);
k1=m(1); k2=m(2)-k1^2; k3=m(3)-3*k1*m(2)+2*k1^3; b=sym(be);
fprintf('near-normal: zeta %.3g skew %.10g sym %.10g nigform %.10g\n', de*ga, r.Skewness, double(vpa((3*b*k2+b^3*k3)/(k1+b^2*k2)^(3/2))), 3*be/(al*sqrt(de*ga)));
% Integer (lattice) data: mid-quantile Bowley vs underlying continuous
d = nigs(rs,P,3000)*3; Bc = bow(d); Bl = bow(round(d));
r = GeneralizedHyperbolic_MLE(round(d),'NIG','NoProfile');
fprintf('lattice: continuous %.4f  plain quantile on rounded %.4f  mid-quantile on rounded %.4f\n', Bc, Bl, r.SkewnessBowleyData);

function b = bow(y), q = quantile(y,[.25 .5 .75]); b=(q(3)+q(1)-2*q(2))/(q(3)-q(1)); end
function x = nigs(rs,p,n)
  m=p.delta/p.gamma; lam=p.delta^2; y=randn(rs,n,1).^2;
  v=m+(m^2*y)/(2*lam)-(m/(2*lam))*sqrt(4*m*lam*y+m^2*y.^2); u=rand(rs,n,1); fl=u>m./(m+v); v(fl)=m^2./v(fl);
  x=p.mu+p.beta*v+sqrt(v).*randn(rs,n,1);
end
function l = ghlog(x,mu,lam,al,be,de)
  ga=sqrt(al^2-be^2); q=sqrt(de^2+(x-mu).^2);
  l = lam*log(ga/de)-0.5*log(2*pi)-(lam-0.5)*log(al)-log(besselk(lam,de*ga)) + (lam-0.5)*log(q) + log(besselk(lam-0.5,al*q,1))-al*q + be*(x-mu);
end
