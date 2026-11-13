scr = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad';
load(fullfile(scr, 'ar1_run', 'ar1_results.mat'), 'allRes');
sets = ["per0dd_sleep_all","per0dd_wake_all","csld_sleep_light","csld_sleep_dark","csld_wake_light","csld_wake_dark"];
xm = [300 2 300 300 2 2];
rs = RandStream('mrg32k3a', 'Seed', 1);
pPrim = NaN(1, 2);
for k = 1:numel(sets)
    A = allRes.(sets(k)); T = A.T; ok = find(isfinite(T.B));
    % formula check against BoutAR1Noise on the real pairs
    Bc = NaN(numel(ok), 1); Hc = Bc;
    for j = 1:numel(ok)
        P = A.R.Residuals{ok(j)};
        [Bc(j), Hc(j)] = bh(P(:, 1), P(:, 2));
    end
    assert(max(abs(Bc - T.B(ok))) < 1e-12 && max(abs(Hc - T.Hetero(ok))) < 1e-12, 'formula mismatch %s', sets(k));
    obs = spear(T.B(ok), T.Hetero(ok));
    n = T.NPairs(ok); mb = T.MeanBout(ok);
    nul = NaN(1000, 1);
    for r = 1:1000
        b = NaN(numel(ok), 1); h = b;
        for j = 1:numel(ok)
            x = xm(k) + round(-(mb(j) - xm(k)) * log(rand(rs, n(j) + 1, 1)));
            [b(j), h(j)] = bh(x(1:end-1), x(2:end));
        end
        nul(r) = spear(b, h);
    end
    p = (1 + sum(nul >= obs)) / 1001;
    % Cox magnitude: latent AR(1) phi 0.5 on log mean, sigma matched to the median B
    medB = median(T.B(ok)); coxR = NaN; sig = NaN;
    if medB > 0
        sg = 0:0.05:2.5; mB = NaN(size(sg));
        for a = 1:numel(sg)
            rs2 = RandStream('mrg32k3a', 'Seed', 2); bb = NaN(numel(ok), 1);
            for j = 1:numel(ok), x = cox(n(j) + 1, mb(j), xm(k), sg(a), rs2); bb(j) = bh(x(1:end-1), x(2:end)); end
            mB(a) = median(bb);
        end
        [~, ia] = min(abs(mB - medB)); sig = sg(ia);
        rc = NaN(200, 1); rs3 = RandStream('mrg32k3a', 'Seed', 3);
        for r = 1:200
            b = NaN(numel(ok), 1); h = b;
            for j = 1:numel(ok), x = cox(n(j) + 1, mb(j), xm(k), sig, rs3); [b(j), h(j)] = bh(x(1:end-1), x(2:end)); end
            rc(r) = spear(b, h);
        end
        coxR = median(rc);
    end
    fprintf('%-18s n %2d  Spearman(B, Hetero) %+.3f  null median %+.3f (95th %.3f)  one-sided p %.4f  | Cox (phi .5, sigma %.2f, median B %.3f): expected %+.3f\n', ...
        sets(k), numel(ok), obs, median(nul), prctile(nul, 95), p, sig, medB, coxR);
    if k <= 2, pPrim(k) = p; end
end
[ps, ix] = sort(pPrim); h = [min(1, 2*ps(1)), max(min(1, 2*ps(1)), ps(2))]; hp(ix) = h;
fprintf('PRIMARY B per0: p sleep %.4f wake %.4f; Holm %.4f %.4f\n', pPrim, hp);

function [B, H] = bh(a, c)
ac = a - mean(a); B = sum(ac .* (c - mean(c))) / sum(ac.^2);
e = c - (mean(c) - B * mean(a)) - B * a;
H = spear(abs(e), a);
end
function r = spear(a, c)
ra = tr(a); rc = tr(c); ra = ra - mean(ra); rc = rc - mean(rc);
r = sum(ra .* rc) / sqrt(sum(ra.^2) * sum(rc.^2));
end
function r = tr(v)
[~, ~, j] = unique(v); cnt = accumarray(j, 1); cum = cumsum(cnt); avg = cum - (cnt - 1) / 2; r = avg(j);
end
function x = cox(n, mb, xmin, sig, rs)
z = zeros(n, 1); z(1) = randn(rs);
for t = 2:n, z(t) = 0.5 * z(t-1) + sqrt(0.75) * randn(rs); end
m = (mb - xmin) * exp(sig * z - sig^2 / 2);
x = xmin + round(-m .* log(rand(rs, n, 1)));
end
