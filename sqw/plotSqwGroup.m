function fh = plotSqwGroup(G)
%PLOTSQWGROUP  Per-fly observable quantities of one or more genotypes (SQWOBSERVABLEGROUP outputs).
%
%   fh = plotSqwGroup(G)
%
%   One panel per quantity (c_obs, count J, count EPR, E3 - E1); one column
%   of jittered points per genotype x phase, in the order given (pass G in
%   the paper's genotype order: CS, per0LD, per0DD, fumin, sss), with the
%   median as a black bar and its fly-bootstrap 95% CI as a grey box. Only
%   fly-phases with a Q/W reference; E3 - E1 only for chains with a cycle.

arguments
    G struct
end
fh = figure('Name', 'sqw observable', 'Color', 'w', 'Position', [100 100 1100 700]);
try, theme(fh, 'light'); catch, end
qs = ["c_obs" "J" "EPR" "E31"];
lab = ["c_{obs}" "count J (cycles/h)" "count EPR (nats/h)" "E_3 - E_1 (nats)"];
for k = 1:4
    ax = subplot(2, 2, k, 'Parent', fh); hold(ax, 'on');
    x = 0; ticks = []; names = strings(0);
    for g = 1:numel(G)
        T = G(g).perFly; S = G(g).summary;
        for p = unique(T.phase)'
            x = x + 1;
            R = T(T.phase == p & T.refOK, :);
            if qs(k) == "E31", R = R(R.hasCycle, :); end
            v = R.(qs(k)); v = v(isfinite(v));
            s = S(S.phase == p, :);
            lo = s.(qs(k) + "_lo"); hi = s.(qs(k) + "_hi"); md = s.(qs(k) + "_med");
            if isfinite(lo), rectangle(ax, 'Position', [x - 0.3, lo, 0.6, max(hi - lo, eps)], ...
                    'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none'); end
            jit = (mod((1:numel(v))' * 0.618, 1) - 0.5) * 0.4;
            scatter(ax, x + jit, v, 14, 'filled', 'MarkerFaceAlpha', 0.6);
            plot(ax, x + [-0.3 0.3], [md md], 'k-', 'LineWidth', 2);
            ticks(end + 1) = x; names(end + 1) = G(g).genotype + " " + s.phaseName; %#ok<AGROW>
        end
    end
    yline(ax, 0, ':');
    xticks(ax, ticks); xticklabels(ax, names); xtickangle(ax, 45);
    ylabel(ax, lab(k)); xlim(ax, [0.4, x + 0.6]); grid(ax, 'on'); box(ax, 'on');
end
end
