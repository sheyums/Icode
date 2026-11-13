icode = 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode';
scr   = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad';
cd(icode); addpath(fileparts(icode), '-end');
for fn = {'process_LocomToSleep','extractBoutsForModelCompetition'}
    fprintf('%s -> %s\n', fn{1}, which(fn{1}));
end
% ---- build tim0 exactly as per0 DD in run_homeo.m
S = process_LocomToSleep(timeInterval=1, segmentMode=1, startZT=0, Teliminate=0, ...
    dayLength=24, daytime=12, nighttime=12, deadThresh=12, minDays=3, bridgeWakeSeconds=1, ...
    folderPath=fullfile(scr,'tim0dd_stage'), selectAll=true);
[Bs, ms] = extractBoutsForModelCompetition(S, 'entire', 'sleep', Verbose=false);
[Bw, mw] = extractBoutsForModelCompetition(S, 'entire', 'wake', Verbose=false);
assert(isequal(ms.DataStartZT, 0) && isequal(mw.DataStartZT, 0), 'not starting at ZT0');
nm = string(S.includedFiles(:));
F = struct('Dur', {}, 'IsSleep', {}, 'Onset', {}, 'Strata', {}, 'Name', {}); drop = strings(0);
for i = 1:size(Bs, 1)
    vs = ~isnan(Bs(i, :)); vw = ~isnan(Bw(i, :));
    if sum(vs) < 40 || sum(vw) < 40, drop(end+1) = nm(i); continue, end %#ok<SAGROW>
    dur = [Bs(i, vs) Bw(i, vw)].'; on = [ms.OnsetBins(i, vs) mw.OnsetBins(i, vw)].' - 1;
    isS = [true(sum(vs), 1); false(sum(vw), 1)];
    [on, o] = sort(on); dur = dur(o); isS = isS(o);
    F(end+1) = struct('Dur', dur, 'IsSleep', isS, 'Onset', on, 'Strata', floor(mod(on / 3600, 24)), 'Name', nm(i)); %#ok<SAGROW>
end
fprintf('tim0 DD: %d files included by pipeline, %d kept (>=40 bouts each type), dropped: %s\n', numel(nm), numel(F), strjoin(drop, ', '));
fprintf('record length h: %s\n', mat2str(unique(round(arrayfun(@(f) (f.Onset(end) + f.Dur(end)) / 3600, F)))));
tim0 = F; save(fullfile(scr, 'N_run', 'tim0_flies.mat'), 'tim0');
L = load(fullfile(scr, 'homeo_run', 'homeo_results.mat'), 'flies'); per0 = L.flies.per0DD;
% ---- statistics
G = {tim0, per0}; gn = ["tim0 DD", "per0 DD (days 1-4, reference)"]; out = struct();
for k = 1:2
    [Lv, Av] = dayhalf(G{k});               % fly x day(4) x half(2)
    dL = Lv(:,1,:); dL = squeeze(mean(Lv(:,1,:),3) - mean(mean(Lv(:,2:4,:),3),2));
    dA = squeeze(sum(Av(:,1,:),3) - mean(sum(Av(:,2:4,:),3),2));
    cL = Lv(:,:,2) - Lv(:,:,1); cA = Av(:,:,2) - Av(:,:,1);
    dcL = cL(:,1) - mean(cL(:,2:4),2); dcA = cA(:,1) - mean(cA(:,2:4),2);
    V = {dL, dA, dcL, dcA}; vn = ["dL level, log bout length", "dA level, sleep min/day", "dcL half-day contrast, log length", "dcA half-day contrast, sleep min"];
    p = cellfun(@(v) flipp(v, 10000, 1), V);
    fprintf('\n#### %s: %d flies\n', gn(k), size(Lv, 1));
    for j = 1:4
        v = V{j}; v = v(isfinite(v));
        fprintf('  %-34s mean %+.3f  (%d/%d flies > 0)  p %.4f\n', vn(j), mean(v), sum(v > 0), numel(v), p(j));
    end
    fprintf('  day-1 contrast alone: cL %+.3f p %.4f; cA %+.1f p %.4f\n', mean(cL(:,1),'omitnan'), flipp(cL(:,1),10000,1), mean(cA(:,1),'omitnan'), flipp(cA(:,1),10000,1));
    fprintf('  profiles (mean over flies), L by day x half:\n'); disp(round(squeeze(mean(Lv, 1, 'omitnan')), 3))
    fprintf('  A (sleep min) by day x half:\n'); disp(round(squeeze(mean(Av, 1, 'omitnan')), 1))
    out(k).name = gn(k); out(k).L = Lv; out(k).A = Av; out(k).V = V; out(k).p = p; %#ok<SAGROW>
end
h = holm(out(1).p);
fprintf('\nPRIMARY tim0 Holm: dL %.4f, dA %.4f, dcL %.4f, dcA %.4f -> DROP first 24 h: %d\n', h, any(h < 0.05));
save(fullfile(scr, 'N_run', 'N_results.mat'), 'out', 'h');

function [Lv, Av] = dayhalf(F)
nF = numel(F); Lv = NaN(nF, 4, 2); Av = NaN(nF, 4, 2);
for f = 1:nF
    o = F(f).Onset(F(f).IsSleep); d = F(f).Dur(F(f).IsSleep);
    for dd = 1:4
        for hh = 1:2
            a = (dd - 1) * 86400 + (hh - 1) * 43200; b = a + 43200;
            k = o >= a & o < b; if sum(k) >= 3, Lv(f, dd, hh) = mean(log(d(k))); end
            Av(f, dd, hh) = sum(max(0, min(o + d, b) - max(o, a))) / 60;
        end
    end
end
end
function p = flipp(v, nR, seed)
v = v(isfinite(v)); rs = RandStream('mrg32k3a', 'Seed', seed);
s = 2 * (rand(rs, numel(v), nR) > 0.5) - 1; n = mean(v .* s, 1); o = mean(v);
p = (1 + sum(abs(n) >= abs(o))) / (nR + 1);
end
function a = holm(p)
[ps, ix] = sort(p); m = numel(p); a = zeros(size(p)); r = 0;
for k = 1:m, r = max(r, min(1, (m - k + 1) * ps(k))); a(ix(k)) = r; end
end
