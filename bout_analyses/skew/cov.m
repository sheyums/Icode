function cov(model, n, R, seed, lam)
% Coverage of SkewnessCIWald; bias of Bowley model/data vs truth.
addpath('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode');
P = struct('mu',0,'lambda',lam,'alpha',2,'beta',0.8,'delta',1); P.gamma = sqrt(P.alpha^2-P.beta^2);
g0 = mskew(P); b0 = bowtrue(P);
rs = RandStream('mt19937ar','Seed',seed);
cv = nan(R,1); sk = cv; se = cv; bm = cv; bd = cv; t = tic;
warning('off','all');
for i = 1:R
  x = ghs(rs,P,n);
  r = GeneralizedHyperbolic_MLE(x, model, 'NoProfile');
  sk(i)=r.Skewness; se(i)=r.SkewnessSE; bm(i)=r.SkewnessBowleyModel; bd(i)=r.SkewnessBowleyData;
  cv(i) = r.SkewnessCIWald(1) <= g0 && g0 <= r.SkewnessCIWald(2);
end
ok = ~isnan(se);
fprintf('%s n=%d lam=%g R=%d | true skew %.4f: coverage %.3f (%d valid) median est %.4f, sd est %.4f, median SE %.4f | Bowley true %.4f: model mean %.4f sd %.4f; data mean %.4f sd %.4f | %.1fs/fit\n', ...
  model,n,lam,R,g0,mean(cv(ok)),nnz(ok),median(sk),std(sk(ok)),median(se(ok)),b0,mean(bm,'omitnan'),std(bm,'omitnan'),mean(bd),std(bd),toc(t)/R);
end
function g = mskew(p)
 w=p.delta*p.gamma; s=p.delta/p.gamma; K=@(v) besselk(v,w);
 m1=s*K(p.lambda+1)/K(p.lambda); m2=s^2*K(p.lambda+2)/K(p.lambda); m3=s^3*K(p.lambda+3)/K(p.lambda);
 k2=m2-m1^2; k3=m3-3*m1*m2+2*m1^3; g=(3*p.beta*k2+p.beta^3*k3)/(m1+p.beta^2*k2)^1.5;
end
function b = bowtrue(p)
 xs=linspace(-80,80,4e6); q=sqrt(p.delta^2+xs.^2);
 l=p.lambda*log(p.gamma/p.delta)-(p.lambda-0.5)*log(p.alpha)-log(besselk(p.lambda,p.delta*p.gamma))+(p.lambda-0.5)*log(q)+log(besselk(p.lambda-0.5,p.alpha*q,1))-p.alpha*q+p.beta*xs;
 F=cumtrapz(xs,exp(l)); F=F/F(end); [Fu,iu]=unique(F); qq=interp1(Fu,xs(iu),[.25 .5 .75]); b=(qq(3)+qq(1)-2*qq(2))/(qq(3)-qq(1));
end
function x = ghs(rs,p,n)
 lf=@(u) p.lambda*u-0.5*(p.delta^2*exp(-u)+p.gamma^2*exp(u));
 u=linspace(-60,60,200001); L=lf(u); m=max(L); k=find(L>m-60);
 u=linspace(u(max(k(1)-1,1)),u(min(k(end)+1,numel(u))),400001);
 F=cumtrapz(u,exp(lf(u)-m)); F=F/F(end); [Fu,iu]=unique(F);
 W=exp(interp1(Fu,u(iu),rand(rs,n,1),'pchip')); x=p.mu+p.beta*W+sqrt(W).*randn(rs,n,1);
end
