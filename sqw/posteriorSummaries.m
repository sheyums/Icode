function S = posteriorSummaries(Z, ph, nP, hours, options)
%POSTERIORSUMMARIES  Per-path statistics averaged over posterior paths, per phase (sqw_ref.posterior_summaries).
%
%   S = posteriorSummaries(Z, ph, nP, hours)
%   S = posteriorSummaries(Z, ph, nP, hours, MinRoutes=3, MinDefined=0.5)
%
%   Z      (nPaths x T) HMMFFBS paths;  hours (nP x 1) hours per phase
%   S      struct, (nP x 1) arrays:
%            Xn, Yn, c   computed on EACH path (XYCOUNTS, >= MinRoutes routes
%                        per denominator), then averaged over the paths where
%                        defined; NaN if defined on < MinDefined of the paths
%            J           mean JCYCLECOUNT (linear, so = J of the mean counts)
%            EPR_bound   mean per-path count EPR (rule-of-three bound)
%            frac_substituted  share of paths where the rule of three acted
%          plus routeCounts (nP x 4) and pairCounts (3 x 3 x nP), means over paths
%   Jensen: nonlinear statistics are averaged per path, E[f(N)], never computed
%   on the mean counts.

arguments
    Z double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    hours (:,1) double {mustBePositive}
    options.MinRoutes (1,1) double = 3
    options.MinDefined (1,1) double = 0.5
end
nZ = size(Z, 1);
X = NaN(nZ, nP); Y = X; J = X; E = X; sub = false(nZ, nP);
rc = zeros(nP, 4, nZ); pc = zeros(3, 3, nP, nZ);
for k = 1:nZ
    z = Z(k, :)';
    C = routeCounts(z, ph, nP); N = pairCounts(z, ph, nP);
    [X(k, :), Y(k, :)] = xyCounts(C, MinRoutes=options.MinRoutes);
    J(k, :) = jCycleCount(N, hours)';
    for p = 1:nP
        e = epr_count_based(N(:, :, p), hours(p) * 3600, 'bound');
        if isnan(e), e = 0; end                        % no edge at all: sqw_ref returns 0
        E(k, p) = e;
        Np = N(:, :, p); u = [Np(1,2) Np(1,3) Np(2,3)]; v = [Np(2,1) Np(3,1) Np(3,2)];
        sub(k, p) = any(xor(u == 0, v == 0));
    end
    rc(:, :, k) = C; pc(:, :, :, k) = N;
end
C = (Y - X) / 2;
S.Xn = pathMean(X, options.MinDefined);
S.Yn = pathMean(Y, options.MinDefined);
S.c = pathMean(C, options.MinDefined);
S.J = mean(J, 1)';
S.EPR_bound = mean(E, 1)';
S.frac_substituted = mean(sub, 1)';
S.routeCounts = mean(rc, 3);
S.pairCounts = mean(pc, 4);
end

function m = pathMean(V, minDef)
ok = mean(isfinite(V), 1) >= minDef;
V(~isfinite(V)) = NaN;
m = mean(V, 1, 'omitnan')';
m(~ok) = NaN;
end
