D = 'C:\Users\sheyu\Box\Data from Shreenath\memory_vs_no_memory_residuals_raw_export';
out = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\F_run';
vars = {'fly_id','replicate','version','residual_seconds'};
pAll = []; lab = strings(0); res = struct();
for ph = ["day","night"]
    tic
    o = rd(fullfile(D, ph, 'residuals.csv'), vars);
    s = rd(fullfile(D, ph, 'sim_residuals.csv'), vars);
    p = rd(fullfile(D, ph, 'sim_residuals.csv'), {'fly_id','generating_phi'});
    fprintf('%s: read %d observed and %d simulated rows in %.0f s\n', ph, height(o), height(s), toc);
    flies = unique(o.fly_id); nF = numel(flies);
    T = table(flies, 'VariableNames', {'fly'});
    T.n = zeros(nF, 1); T.phi = NaN(nF, 1);
    for k = ["obs_","fit_","p0_"], for st = ["skew","bowley","sd"], T.(k + st) = NaN(nF, 1); end, end
    T.p_adequacy = NaN(nF, 1); T.fit_skew_sd = NaN(nF, 1);
    [gs, gf, gv, gr] = findgroups(s.fly_id, s.version, s.replicate);
    S3 = splitapply(@(x) {stats3(x)}, s.residual_seconds, gs);
    S3 = vertcat(S3{:});
    for f = 1:nF
        x = o.residual_seconds(o.fly_id == flies(f));
        T.n(f) = numel(x); so = stats3(x);
        T.obs_skew(f) = so(1); T.obs_bowley(f) = so(2); T.obs_sd(f) = so(3);
        T.phi(f) = median(p.generating_phi(p.fly_id == flies(f)), 'omitnan');
        mf = gf == flies(f) & gv == "fitted"; m0 = gf == flies(f) & gv == "phi0";
        T.fit_skew(f) = mean(S3(mf, 1), 'omitnan'); T.fit_bowley(f) = mean(S3(mf, 2), 'omitnan'); T.fit_sd(f) = mean(S3(mf, 3), 'omitnan');
        T.p0_skew(f) = mean(S3(m0, 1), 'omitnan'); T.p0_bowley(f) = mean(S3(m0, 2), 'omitnan'); T.p0_sd(f) = mean(S3(m0, 3), 'omitnan');
        nul = S3(mf, 1); nul = nul(isfinite(nul)); T.fit_skew_sd(f) = std(nul);
        T.p_adequacy(f) = min(1, 2 * min(1 + sum(nul >= so(1)), 1 + sum(nul <= so(1))) / (numel(nul) + 1));
    end
    nsig = sum(T.p_adequacy < 0.05); pc = binom_upper(nsig, nF, 0.05);
    d = T.fit_skew - T.p0_skew; rs = RandStream('mrg32k3a', 'Seed', 1);
    fl = mean(d .* (2 * (rand(rs, nF, 10000) > 0.5) - 1), 1);
    pm = min(1, 2 * min(1 + sum(fl >= mean(d)), 1 + sum(fl <= mean(d))) / 10001);
    fprintf(['\n=== %s, %d flies; phi median %.3f (IQR %.3f to %.3f)\n' ...
        '  (1) adequacy: %d/%d flies outside their fitted-model skewness, count p %.4f\n' ...
        '  (2) memory: mean skew fitted %.3f vs phi0 %.3f, diff %.4f (%.1f%% of phi0), sign-flip p %.4f\n' ...
        '  observed skew mean %.3f; Bowley obs/fit/phi0 %.3f/%.3f/%.3f; SD obs/fit/phi0 %.0f/%.0f/%.0f s\n'], ...
        ph, nF, median(T.phi), prctile(T.phi, 25), prctile(T.phi, 75), nsig, nF, pc, mean(T.fit_skew), mean(T.p0_skew), ...
        mean(d), 100 * mean(d) / mean(T.p0_skew), pm, mean(T.obs_skew), mean(T.obs_bowley), mean(T.fit_bowley), ...
        mean(T.p0_bowley), mean(T.obs_sd), mean(T.fit_sd), mean(T.p0_sd));
    db = T.fit_bowley - T.p0_bowley; ds = T.fit_sd - T.p0_sd;
    fprintf('  secondary: Bowley fit-phi0 mean %.4f; SD fit-phi0 mean %.1f s; flies with fit skew < phi0 skew: %d/%d\n', mean(db), mean(ds), sum(d < 0), nF);
    pAll = [pAll pc pm]; lab = [lab ph + " adequacy", ph + " memory"]; %#ok<AGROW>
    res.(ph) = T; writetable(T, fullfile(out, "F_" + ph + "_perfly.csv"));
end
[ps, ix] = sort(pAll); m = numel(pAll); h = zeros(1, m); run = 0;
for k = 1:m, run = max(run, min(1, (m - k + 1) * ps(k))); h(ix(k)) = run; end
fprintf('\nPRIMARY (Holm across 4):\n'); for k = 1:m, fprintf('  %-16s p %.4f  Holm %.4f\n', lab(k), pAll(k), h(k)); end
save(fullfile(out, 'F_results.mat'), 'res');

function v = stats3(x)
x = x(isfinite(x)); n = numel(x); v = NaN(1, 3);
if n < 10, return, end
xc = x - mean(x); m2 = mean(xc.^2);
v(1) = mean(xc.^3) / m2^1.5; v(3) = sqrt(m2);
[u, ~, j] = unique(sort(x)); c = accumarray(j, 1); Fm = (cumsum(c) - c / 2) / n;
P = [0.25 0.5 0.75];
if numel(u) > 1
    q = interp1(Fm, u, P, 'linear'); q(P < Fm(1)) = u(1); q(P > Fm(end)) = u(end);
    if q(3) > q(1), v(2) = (q(3) + q(1) - 2 * q(2)) / (q(3) - q(1)); end
end
end
function p = binom_upper(k, n, q)
j = k:n; lp = gammaln(n + 1) - gammaln(j + 1) - gammaln(n - j + 1) + j * log(q) + (n - j) * log1p(-q);
p = min(1, sum(exp(lp)));
end
function t = rd(f, vars)
op = detectImportOptions(f, 'TextType', 'string');
op.SelectedVariableNames = vars;
t = readtable(f, op);
end
