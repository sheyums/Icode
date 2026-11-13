icode = 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode';
scr   = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad';
out   = fullfile(scr, 'homeo_run');
cd(icode); addpath(fileparts(icode), '-end');
for fn = {'BoutCrossType','SleepAmountACF','ProcessSProxy','RenewalSurrogate'}
    assert(isequal(fileparts(which(fn{1})), pwd), '%s', fn{1});
end
% ---- datasets: {name, {folder, startZT, Teliminate} ..., profile?}
D = struct('name', {'per0DD','CSLD','CSDD','per0LD'}, ...
    'parts', {{{fullfile(scr,'per0dd_stage'),0,0}}, {{fullfile(scr,'csld_stage'),0,0}}, ...
              {{fullfile(scr,'csdd_stage'),0,0}}, {{fullfile(scr,'per0ld_stage','8_25_2017'),20,4}}}, ...
    'profile', {false, true, true, true});
flies = struct();
for d = D
    F = struct('Dur', {}, 'IsSleep', {}, 'Onset', {}, 'Strata', {}, 'Name', {});
    for p = 1:numel(d.parts)
        pt = d.parts{p};
        S = process_LocomToSleep(timeInterval=1, segmentMode=1, startZT=pt{2}, Teliminate=pt{3}, ...
            dayLength=24, daytime=12, nighttime=12, deadThresh=12, minDays=3, bridgeWakeSeconds=1, ...
            folderPath=pt{1}, selectAll=true);
        [Bs, ms] = extractBoutsForModelCompetition(S, 'entire', 'sleep', Verbose=false);
        [Bw, mw] = extractBoutsForModelCompetition(S, 'entire', 'wake', Verbose=false);
        assert(isequal(ms.DataStartZT, 0) && isequal(mw.DataStartZT, 0), 'not starting at ZT0');
        nm = string(S.includedFiles(:));
        for i = 1:size(Bs, 1)
            vs = ~isnan(Bs(i, :)); vw = ~isnan(Bw(i, :));
            if sum(vs) < 40 || sum(vw) < 40, continue, end
            dur = [Bs(i, vs) Bw(i, vw)].'; on = [ms.OnsetBins(i, vs) mw.OnsetBins(i, vw)].' - 1;
            isS = [true(sum(vs), 1); false(sum(vw), 1)];
            [on, o] = sort(on); dur = dur(o); isS = isS(o);
            F(end+1) = struct('Dur', dur, 'IsSleep', isS, 'Onset', on, ...
                'Strata', floor(mod(on / 3600, 24)), 'Name', nm(i)); %#ok<SAGROW>
        end
    end
    nadj = arrayfun(@(f) sum(f.Onset(2:end) ~= f.Onset(1:end-1) + f.Dur(1:end-1)), F);
    alt = arrayfun(@(f) all(f.IsSleep(2:end) ~= f.IsSleep(1:end-1)), F);
    fprintf('\n#### %s: %d flies; non-adjacent links per fly: median %g, max %d; strictly alternating: %d/%d\n', ...
        d.name, numel(F), median(nadj), max(nadj), sum(alt), numel(F));
    flies.(d.name) = F;
end
% ---- A
res = struct();
pA = [];
for d = D
    F = flies.(d.name);
    if d.profile, nul = "strata"; else, nul = "shift"; end
    R = BoutCrossType(F, Null=nul, NumSurrogates=1000, Seed=1, MinTriples=20, FlyNames=[F.Name]);
    res.A.(d.name) = R;
    fprintf('\n=== A %s (%s null), %d flies\n', d.name, nul, R.Group.NFlies(1)); disp(R.Group)
end
G = res.A.per0DD.Group;
p = [G.MeanP(G.Stat == "r_pre"), G.MeanP(G.Stat == "Delta")];
fprintf('PRIMARY A per0 DD: r_pre mean %.4f p %.4f; Delta mean %.4f p %.4f; Holm %.4f %.4f\n', ...
    G.Mean(G.Stat == "r_pre"), p(1), G.Mean(G.Stat == "Delta"), p(2), holm(p));
pc = arrayfun(@(n) res.A.(n).Group.CountP(res.A.(n).Group.Stat == "r_pre"), ["CSLD","CSDD","per0LD"]);
nc = arrayfun(@(n) res.A.(n).Group.NSig(res.A.(n).Group.Stat == "r_pre"), ["CSLD","CSDD","per0LD"]);
fprintf('PRIMARY A profile sets r_pre counts CSLD %d, CSDD %d, per0LD %d; count p %s; Holm %s\n', nc, mat2str(pc, 4), mat2str(holm(pc), 4));
% ---- C (descriptive)
for bm = [30 60]
    R = SleepAmountACF(flies.per0DD, BinMinutes=bm, NumSurrogates=1000, Seed=1);
    res.C.("bin" + bm) = R;
    fprintf('\n=== C per0 DD, %d-min windows (DESCRIPTIVE)\n', bm); disp(R.Group)
end
% ---- D
R = ProcessSProxy(flies.per0DD, NumSurrogates=1000, Seed=1);
res.D = R;
fprintf('\n=== D per0 DD\n'); disp(R.Group)
fprintf('group rho_S grid (rows TauW 0.5..16, cols TauS 0.5..16):\n'); disp(round(R.GroupRhoS, 3))
fprintf('group rho_W grid:\n'); disp(round(R.GroupRhoW, 3))
rA = res.A.per0DD.PerFly.r_pre; rD = R.PerFly.RhoS_atGroupBest;
fprintf('per fly: Spearman(rho_S at group best, r_pre) = %.3f over %d flies\n', spear(rA, rD), sum(isfinite(rA) & isfinite(rD)));
save(fullfile(out, 'homeo_results.mat'), 'res', 'flies', '-v7.3');
fprintf('DONE\n');

function a = holm(p)
[ps, ix] = sort(p); m = numel(p); a = zeros(size(p)); run = 0;
for k = 1:m, run = max(run, min(1, (m - k + 1) * ps(k))); a(ix(k)) = run; end
end
function r = spear(a, b)
ok = isfinite(a) & isfinite(b); a = tiedrank(a(ok)); b = tiedrank(b(ok));
r = corr(a, b);
end
