icode  = 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode';
parent = fileparts(icode);
stage  = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\csld_stage';
out    = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\csld_run';
cd(icode); addpath(parent, '-end');
assert(isequal(fileparts(which('BoutIrreversibility')), pwd));
fprintf('process_LocomToSleep: %s\n', which('process_LocomToSleep'));
assert(numel(dir(fullfile(stage, '*.txt'))) == 40);
S = process_LocomToSleep(timeInterval=1, segmentMode=1, Teliminate=0, startZT=0, ...
    dayLength=24, daytime=12, nighttime=12, deadThresh=12, minDays=3, bridgeWakeSeconds=1, ...
    folderPath=stage, selectAll=true);
fprintf('flies retained by the pipeline: %d of 40\n', numel(S.includedFiles));
types = ["sleep","wake"];
for ti = 1:2
    [B, meta] = extractBoutsForModelCompetition(S, 'entire', char(types(ti)), Verbose=false);
    fprintf('DataStartZT = %s\n', mat2str(meta.DataStartZT));
    assert(isequal(meta.DataStartZT, 0), 'records do not start at ZT0');
    O = meta.OnsetBins;
    for i = 1:size(O, 1)
        o = O(i, ~isnan(O(i, :)));
        assert(all(diff(o) > 0), 'bouts of fly %d not in time order', i);
    end
    ZT = mod((O - 1) / 3600, 24);
    St = floor(ZT);
    nb = sum(~isnan(B), 2);
    keep = nb >= 40;
    names = string(S.includedFiles(:));
    fprintf('\n===== %s: %d flies in, %d with >= 40 bouts; bouts per kept fly %d-%d (median %g); min bout %g s\n', ...
        types(ti), numel(nb), sum(keep), min(nb(keep)), max(nb(keep)), median(nb(keep)), meta.MinObservedSeconds);
    if any(~keep), fprintf('  dropped: %s\n', strjoin(names(~keep) + " (" + nb(~keep) + ")", ', ')); end
    tic
    R = BoutIrreversibility(B(keep, :), Strata=St(keep, :), Lags=[1 2 4 8], NumSurrogates=1000, ...
        Null="both", Seed=1, MinPairs=20, FlyNames=names(keep).', ReturnNull=true);
    fprintf('  BoutIrreversibility: %.0f s\n', toc);
    G = R.Group;
    disp(G(:, {'Lag','NFlies','UpFraction','UpFraction_NSig_Strata','UpFraction_CountP_Strata', ...
        'UpFraction_P_IAAFT','UpFraction_P_Shuffle','Bowley','Rho'}))
    T = R.PerFly;
    for k = [1 2 4 8]
        r = T(T.Lag == k, :);
        nullUp = squeeze(mean(R.PerFlyNull.strata.UpFraction(:, R.Settings.Lags == k, :), 3, 'omitnan'));
        sig = r.UpFraction_P_Strata < 0.05;
        fprintf(['  lag %d: Up strata-sig %d (above own null %d, below %d); Bowley strata-sig %d; ' ...
            'Rho strata-sig %d (positive %d); movable %.2f-%.2f; Up minus strata-null mean: median %+.4f\n'], ...
            k, sum(sig), sum(sig & r.UpFraction > nullUp), sum(sig & r.UpFraction < nullUp), ...
            sum(r.Bowley_P_Strata < 0.05), sum(r.Rho_P_Strata < 0.05), sum(r.Rho_P_Strata < 0.05 & r.Rho > 0), ...
            min(r.StrataMovable), max(r.StrataMovable), median(r.UpFraction - nullUp));
    end
    res.(types(ti)) = R; res.(types(ti) + "_bouts") = B(keep, :); res.(types(ti) + "_strata") = St(keep, :);
    res.(types(ti) + "_meta") = meta;
    writetable(T, fullfile(out, "csld_" + types(ti) + "_perfly.csv"));
end
lab = ["sleep","wake"];
p = [res.sleep.Group.UpFraction_CountP_Strata(1), res.wake.Group.UpFraction_CountP_Strata(1)];
adj = min(1, [2 1] .* sort(p)); adj(2) = max(adj); [~, ix] = sort(p); hp(ix) = adj;
fprintf('\n===== PRIMARY (lag 1, UpFraction count vs Binomial(n, 0.05), Holm across 2)\n');
for k = 1:2
    fprintf('  %-6s %d of %d flies significant, count p = %.4f, Holm = %.4f\n', lab(k), ...
        res.(lab(k)).Group.UpFraction_NSig_Strata(1), res.(lab(k)).Group.NFlies(1), p(k), hp(k));
end
res.primary = table(lab.', p.', hp.', 'VariableNames', {'Type','CountP','HolmP'});
save(fullfile(out, 'csld_results.mat'), 'res', '-v7.3');
