L = load('C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\homeo_run\homeo_results.mat', 'flies');
nR = 1000; nulls = ["N1_1h","N2_seg","N3_plain","N1_2h","N1_3h"];
sets = struct('name', {"CSLD","per0DD"}, 'nph', {2, 1});
res = struct(); pPrim = []; labP = strings(0);
for st = sets
    F = L.flies.(st.name); nF = numel(F);
    for p = 1:st.nph
        obs = NaN(nF, 2); nul = NaN(nF, nR, 2, numel(nulls)); same = NaN(nF, 1); ngap = NaN(nF, 1);
        for f = 1:nF
            fl = F(f);
            ngap(f) = sum(fl.Onset(2:end) ~= fl.Onset(1:end-1) + fl.Dur(1:end-1));
            on = fl.Onset(fl.IsSleep); d = fl.Dur(fl.IsSleep);
            zt = mod(on, 86400) / 3600; dy = floor(on / 86400);
            if st.nph == 2, ph = double(zt >= 12) + 1; seg = 2 * dy + ph; else, ph = ones(size(on)); seg = floor(on / 43200); end
            k = ph == p; d = d(k); seg = seg(k); zt = zt(k);
            obs(f, :) = cstat(d, seg);
            pr = seg(2:end) == seg(1:end-1);
            same(f) = mean(floor(zt([false; pr])) == floor(zt([pr; false])));
            G = {floor(zt), seg, ones(size(d)), floor(zt / 2), floor(zt / 3)};
            rs = RandStream('mrg32k3a', 'Seed', 1); rs.Substream = f;
            for j = 1:numel(nulls)
                for r = 1:nR
                    nul(f, r, :, j) = permute(cstat(pw(d, G{j}, rs), seg), [1 3 2]);
                end
            end
        end
        phn = ["day","night"]; pn = st.name + "_" + phn(p);
        if st.nph == 1, pn = st.name; end
        fprintf('\n#### %s: %d flies; median gaps per fly %g; consecutive pairs sharing a 1-h bin same day: %.2f\n', pn, nF, median(ngap), mean(same, 'omitnan'));
        tn = ["upper","lower"];
        for t = 1:2
            go = mean(obs(:, t), 'omitnan');
            fprintf('  %s 10%%: observed %.3f |', tn(t), go);
            for j = 1:numel(nulls)
                g = squeeze(mean(nul(:, :, t, j), 1, 'omitnan'));
                pj = tp(go, g);
                fprintf(' %s %.3f p %.4f |', nulls(j), mean(g), pj);
                if t == 1 && ((st.nph == 2 && j == 1) || (st.nph == 1 && j == 3))
                    pPrim(end+1) = pj; labP(end+1) = pn + " upper vs " + nulls(j); %#ok<AGROW>
                end
            end
            fprintf('\n');
        end
        res.(pn) = struct('obs', obs, 'nul', nul, 'nulls', nulls, 'same', same);
    end
end
[ps, ix] = sort(pPrim); m = numel(pPrim); h = zeros(1, m); r = 0;
for k = 1:m, r = max(r, min(1, (m - k + 1) * ps(k))); h(ix(k)) = r; end
fprintf('\nPRIMARY Holm: '); for k = 1:m, fprintf('%s p %.4f Holm %.4f; ', labP(k), pPrim(k), h(k)); end; fprintf('\n');
save(fullfile(fileparts(mfilename('fullpath')), 'I_results.mat'), 'res', 'pPrim', 'h', 'labP', '-v7.3');

function x = pw(x, g, rs)
[~, ~, gi] = unique(g);
for u = 1:max(gi)
    i = find(gi == u);
    if numel(i) > 1, x(i) = x(i(randperm(rs, numel(i)))); end
end
end
function c = cstat(d, seg)
c = NaN(1, 2); n = numel(d); if n < 10, return, end
sd = sort(d); pr = find(seg(2:end) == seg(1:end-1));
e1 = d >= sd(ceil(0.9 * n)); e2 = d <= sd(ceil(0.1 * n)); E = {e1, e2};
for k = 1:2
    e = E{k}; s = pr(e(pr)); if isempty(s), continue, end
    c(k) = mean(e(s + 1)) / mean(e);
end
end
function p = tp(o, n)
n = n(isfinite(n)); p = min(1, 2 * min(1 + sum(n >= o), 1 + sum(n <= o)) / (numel(n) + 1));
end
