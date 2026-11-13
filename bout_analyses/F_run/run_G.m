D = 'C:\Users\sheyu\Box\Data from Shreenath\memory_vs_no_memory_residuals_raw_export';
vars = {'fly_id','replicate','version','latent_residual','generating_phi'};
pAll = []; lab = strings(0);
for ph = ["day","night"]
    o = rd(fullfile(D, ph, 'residuals.csv'), vars(1:4));
    s = rd(fullfile(D, ph, 'sim_residuals.csv'), vars);
    flies = unique(o.fly_id); nF = numel(flies);
    [gs, gf, gv] = findgroups(s.fly_id, s.version, s.replicate);
    M = splitapply(@(x) {mom(x)}, s.latent_residual, gs); M = vertcat(M{:});
    gv = string(gv); gf = string(gf);
    phiF = splitapply(@(x) median(x, 'omitnan'), s.generating_phi, findgroups(s.fly_id));
    [~, ~, jf] = unique(s.fly_id); ufl = unique(s.fly_id);
    P = NaN(nF, 2); obs = NaN(nF, 3); fit = obs; p0 = obs; phi = NaN(nF, 1);
    for f = 1:nF
        x = o.latent_residual(o.fly_id == flies(f)); obs(f, :) = mom(x);
        mf = gf == flies(f) & gv == "fitted"; m0 = gf == flies(f) & gv == "phi0";
        fit(f, :) = mean(M(mf, :), 1, 'omitnan'); p0(f, :) = mean(M(m0, :), 1, 'omitnan');
        phi(f) = phiF(ufl == flies(f));
        for k = 1:2
            nl = M(mf, k); nl = nl(isfinite(nl));
            P(f, k) = min(1, 2 * min(1 + sum(nl >= obs(f, k)), 1 + sum(nl <= obs(f, k))) / (numel(nl) + 1));
        end
    end
    for k = 1:2
        nsig = sum(P(:, k) < 0.05); pAll(end+1) = bu(nsig, nF, 0.05); %#ok<AGROW>
        nm = ["skewness","excess kurtosis"]; lab(end+1) = ph + " " + nm(k); %#ok<AGROW>
        fprintf('%s %s: %d/%d flies outside, count p %.4f\n', ph, nm(k), nsig, nF, pAll(end));
    end
    fprintf(['  %s means obs/fitted/phi0: skew %.3f/%.3f/%.3f; excess kurtosis %.3f/%.3f/%.3f; SD %.3f/%.3f/%.3f;' ...
        ' sqrt(1-phi^2) mean %.4f\n'], ph, mean(obs(:,1)), mean(fit(:,1)), mean(p0(:,1)), mean(obs(:,2)), ...
        mean(fit(:,2)), mean(p0(:,2)), mean(obs(:,3)), mean(fit(:,3)), mean(p0(:,3)), mean(sqrt(1 - phi.^2)));
end
[ps, ix] = sort(pAll); m = numel(pAll); h = zeros(1, m); r = 0;
for k = 1:m, r = max(r, min(1, (m - k + 1) * ps(k))); h(ix(k)) = r; end
fprintf('PRIMARY Holm: '); for k = 1:m, fprintf('%s %.4f; ', lab(k), h(k)); end; fprintf('\n');
function v = mom(x)
x = x(isfinite(x)); v = NaN(1, 3); if numel(x) < 10, return, end
c = x - mean(x); m2 = mean(c.^2); v = [mean(c.^3) / m2^1.5, mean(c.^4) / m2^2 - 3, sqrt(m2)];
end
function p = bu(k, n, q)
j = k:n; p = min(1, sum(exp(gammaln(n+1) - gammaln(j+1) - gammaln(n-j+1) + j*log(q) + (n-j)*log1p(-q))));
end
function t = rd(f, vars)
op = detectImportOptions(f, 'TextType', 'string'); op.SelectedVariableNames = vars; t = readtable(f, op);
end
