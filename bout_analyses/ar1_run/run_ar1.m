icode  = 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode';
parent = fileparts(icode);
scr    = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad';
out    = fullfile(scr, 'ar1_run');
cd(icode); addpath(parent, '-end');
for fn = {'BoutAR1Noise','FitResidualFamilies','GeneralizedHyperbolic_MLE','AsymLaplace_MLE', ...
          'shiftlognormal_MLE','shiftweibull_MLE','GenNorm_MLE'}
    assert(isequal(fileparts(which(fn{1})), pwd), '%s does not resolve to Icode', fn{1});
end
ds = struct('name', {'per0dd','csld'}, 'stage', {fullfile(scr,'per0dd_stage'), fullfile(scr,'csld_stage')}, ...
    'nfiles', {48, 40});
xminOf = struct('sleep', 300, 'wake', 2);
allRes = struct();
for d = ds
    assert(numel(dir(fullfile(d.stage, '*.txt'))) == d.nfiles);
    S = process_LocomToSleep(timeInterval=1, segmentMode=1, Teliminate=0, startZT=0, dayLength=24, ...
        daytime=12, nighttime=12, deadThresh=12, minDays=3, bridgeWakeSeconds=1, ...
        folderPath=d.stage, selectAll=true);
    fprintf('\n##### %s: %d of %d flies retained\n', d.name, numel(S.includedFiles), d.nfiles);
    names = string(S.includedFiles(:));
    for type = ["sleep","wake"]
        [B, meta] = extractBoutsForModelCompetition(S, 'entire', char(type), Verbose=false);
        assert(isequal(meta.DataStartZT, 0));
        O = meta.OnsetBins;
        for i = 1:size(O, 1)
            o = O(i, ~isnan(O(i, :))); assert(all(diff(o) > 0));
        end
        keep = find(sum(~isnan(B), 2) >= 40);
        xmin = xminOf.(type);
        if d.name == "per0dd"
            sets = {"all"};
        else
            sets = {"light","dark"};
        end
        for si = 1:numel(sets)
            set = sets{si};
            if set == "all"
                bouts = B(keep, :); strata = [];
            else
                bouts = cell(numel(keep), 1); strata = bouts;
                for k = 1:numel(keep)
                    i = keep(k); v = ~isnan(B(i, :));
                    x = B(i, v).'; on = O(i, v).';
                    zt = mod(on - 1, 86400) / 3600; day = floor((on - 1) / 86400);
                    inPh = (zt < 12) == (set == "light");
                    days = unique(day(inPh)).';
                    bouts{k} = arrayfun(@(dd) x(inPh & day == dd), days, 'UniformOutput', false);
                    strata{k} = arrayfun(@(dd) floor(zt(inPh & day == dd)), days, 'UniformOutput', false);
                end
            end
            tag = d.name + "_" + type + "_" + set;
            tic
            if isempty(strata)
                R = BoutAR1Noise(bouts, Xmin=xmin, NumSurrogates=1000, Seed=1, MinPairs=20, ...
                    ReturnResiduals=true, FlyNames=names(keep).');
            else
                R = BoutAR1Noise(bouts, Xmin=xmin, Strata=strata, NumSurrogates=1000, Seed=1, ...
                    MinPairs=20, ReturnResiduals=true, FlyNames=names(keep).');
            end
            t1 = toc;
            T = R.PerFly;
            fam = strings(height(T), 1); dSE = NaN(height(T), 1); famTabs = cell(height(T), 1);
            for f = 1:height(T)
                if isempty(R.Residuals{f}), continue, end
                FT = FitResidualFamilies(R.Residuals{f}(:, 3));
                famTabs{f} = FT;
                fam(f) = FT.Family(find(FT.DeltaBIC == 0, 1));
                dSE(f) = FT.DeltaBIC(FT.Family == "ShiftedExp");
            end
            T.BICWinner = fam; T.ShiftedExpDeltaBIC = dSE;
            fprintf('\n=== %s: %d flies fitted (%.0f s + families %.0f s)\n', tag, sum(isfinite(T.B)), t1, toc - t1);
            disp(R.Group)
            ok = isfinite(T.B);
            [u, ~, j] = unique(fam(ok)); c = accumarray(j, 1);
            fprintf('  BIC winners: %s; ShiftedExp DeltaBIC median %.1f\n', strjoin(u + " " + string(c), ', '), median(dSE(ok)));
            fprintf('  B median %.3f (IQR %.3f to %.3f); Hetero median %.3f; GenBelowXmin median %.4f, max %.4f; GenKS %.3f vs ref %.3f\n', ...
                median(T.B(ok)), prctile(T.B(ok), 25), prctile(T.B(ok), 75), median(T.Hetero(ok)), ...
                median(T.GenBelowXmin(ok)), max(T.GenBelowXmin(ok)), median(T.GenKS(ok)), median(T.GenKS_Ref(ok)));
            fprintf('  Resid skew median %.2f (no-dynamics %.2f); Bowley %.3f (%.3f); SD/mean bout %.3f (%.3f)\n', ...
                median(T.ResidSkew(ok)), median(T.ResidSkew_Null(ok)), median(T.ResidBowley(ok)), ...
                median(T.ResidBowley_Null(ok)), median(T.ResidSD(ok) ./ (T.MeanBout(ok) - xmin)), ...
                median(T.ResidSD_Null(ok) ./ (T.MeanBout(ok) - xmin)));
            writetable(T, fullfile(out, tag + "_perfly.csv"));
            allRes.(tag).R = R; allRes.(tag).T = T; allRes.(tag).families = famTabs;
        end
    end
end
% Holm within dataset on the Hetero count p
for dn = ["per0dd","csld"]
    tags = string(fieldnames(allRes)); tags = tags(startsWith(tags, dn));
    p = arrayfun(@(t) allRes.(t).R.Group.Hetero_CountP, tags);
    [ps, ix] = sort(p); m = numel(p); adj = zeros(m, 1); run = 0;
    for k = 1:m, run = max(run, min(1, (m - k + 1) * ps(k))); adj(ix(k)) = run; end
    fprintf('\n##### PRIMARY %s (Hetero count, Holm across %d)\n', dn, m);
    for k = 1:m
        G = allRes.(tags(k)).R.Group;
        fprintf('  %-22s %2d of %2d flies, count p %.4f, Holm %.4f\n', tags(k), G.Hetero_NSig, G.NFlies, p(k), adj(k));
    end
end
save(fullfile(out, 'ar1_results.mat'), 'allRes', '-v7.3');
fprintf('DONE\n');
