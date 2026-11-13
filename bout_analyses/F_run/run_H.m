D = 'C:\Users\sheyu\Box\Data from Shreenath\memory_vs_no_memory_residuals_raw_export';
vars = {'fly_id','replicate','version','block_index','bout_index_within_block','observed_duration_seconds'};
pAll = []; lab = strings(0); out = struct();
for ph = ["day","night"]
    o = rd(fullfile(D, ph, 'residuals.csv'), vars);
    s = rd(fullfile(D, ph, 'sim_residuals.csv'), vars);
    flies = unique(o.fly_id); nF = numel(flies);
    s = sortrows(s, {'fly_id','version','replicate','block_index','bout_index_within_block'});
    [gs, gf, gv, gr] = findgroups(s.fly_id, s.version, s.replicate);
    C = splitapply(@(b, i, d) {cstat(b, i, d)}, s.block_index, s.bout_index_within_block, ...
        s.observed_duration_seconds, gs); C = vertcat(C{:});          % [U10 L10 U20 L20]
    gf = string(gf); gv = string(gv);
    obs = NaN(nF, 4); fitR = NaN(nF, 200, 4); p0R = fitR; shR = NaN(nF, 1000, 4); P = NaN(nF, 2);
    rs = RandStream('mrg32k3a', 'Seed', 1);
    for f = 1:nF
        x = sortrows(o(o.fly_id == flies(f), :), {'block_index','bout_index_within_block'});
        obs(f, :) = cstat(x.block_index, x.bout_index_within_block, x.observed_duration_seconds);
        mf = gf == flies(f) & gv == "fitted"; m0 = gf == flies(f) & gv == "phi0";
        fitR(f, gr(mf) + 1, :) = permute(C(mf, :), [3 1 2]);
        p0R(f, gr(m0) + 1, :) = permute(C(m0, :), [3 1 2]);
        rs.Substream = f;
        for r = 1:1000
            shR(f, r, :) = permute(cstat(x.block_index, x.bout_index_within_block, ...
                x.observed_duration_seconds(randperm(rs, height(x)))), [3 1 2]);
        end
        for k = 1:2
            nl = fitR(f, :, k); nl = nl(isfinite(nl));
            if isfinite(obs(f, k)), P(f, k) = tp(obs(f, k), nl); end
        end
    end
    tn = ["upper 10%","lower 10%","upper 20%","lower 20%"];
    for k = 1:4
        go = mean(obs(:, k), 'omitnan');
        gfit = squeeze(mean(fitR(:, :, k), 1, 'omitnan'));
        g0 = squeeze(mean(p0R(:, :, k), 1, 'omitnan'));
        gsh = squeeze(mean(shR(:, :, k), 1, 'omitnan'));
        pg = tp(go, gfit);
        fprintf('%s %s: obs %.3f (%d flies) | fitted %.3f [2.5-97.5%%: %.3f-%.3f] p %.4f | phi0 %.3f | shuffle %.3f\n', ...
            ph, tn(k), go, sum(isfinite(obs(:, k))), mean(gfit), prctile_(gfit, 2.5), prctile_(gfit, 97.5), pg, mean(g0), mean(gsh));
        if k <= 2
            pAll(end+1) = pg; lab(end+1) = ph + " " + tn(k); %#ok<AGROW>
            n = sum(isfinite(P(:, k))); ns = sum(P(:, k) < 0.05);
            fprintf('   per fly: %d/%d outside fitted range, Binomial p %.4f; obs>fit mean in %d flies\n', ...
                ns, n, bu(ns, n, 0.05), sum(obs(:, k) > mean(fitR(:, :, k), 2, 'omitnan')));
        end
    end
    out.(ph) = struct('obs', obs, 'fitR', fitR, 'p0R', p0R, 'shR', shR, 'P', P, 'flies', flies);
end
[ps, ix] = sort(pAll); m = numel(pAll); h = zeros(1, m); r = 0;
for k = 1:m, r = max(r, min(1, (m - k + 1) * ps(k))); h(ix(k)) = r; end
fprintf('PRIMARY Holm: '); for k = 1:m, fprintf('%s %.4f; ', lab(k), h(k)); end; fprintf('\n');
save(fullfile(fileparts(mfilename('fullpath')), 'H_results.mat'), 'out', 'pAll', 'h', 'lab');

function c = cstat(b, i, d)
c = NaN(1, 4); n = numel(d); if n < 10, return, end
sd = sort(d);
pr = find(b(2:end) == b(1:end-1) & i(2:end) == i(1:end-1) + 1);   % pair starts
qs = [0.9 0.1 0.8 0.2];
for k = 1:4
    if mod(k, 2) == 1, e = d >= sd(ceil(qs(k) * n)); else, e = d <= sd(ceil(qs(k) * n)); end
    st = pr(e(pr));
    if isempty(st), continue, end
    c(k) = mean(e(st + 1)) / mean(e);
end
end
function p = tp(o, n)
n = n(isfinite(n)); p = min(1, 2 * min(1 + sum(n >= o), 1 + sum(n <= o)) / (numel(n) + 1));
end
function v = prctile_(x, q)
x = sort(x(isfinite(x))); v = x(max(1, round(q / 100 * numel(x))));
end
function p = bu(k, n, q)
j = k:n; p = min(1, sum(exp(gammaln(n+1) - gammaln(j+1) - gammaln(n-j+1) + j*log(q) + (n-j)*log1p(-q))));
end
function t = rd(f, vars)
op = detectImportOptions(f, 'TextType', 'string'); op.SelectedVariableNames = vars; t = readtable(f, op);
end
