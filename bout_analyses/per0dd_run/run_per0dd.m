icode  = 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode';
parent = fileparts(icode);
stage  = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\per0dd_stage';
out    = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\per0dd_run';
cd(icode); addpath(parent, '-end');
assert(isequal(fileparts(which('BoutIrreversibility')), pwd));
fprintf('process_LocomToSleep: %s\n', which('process_LocomToSleep'));
assert(numel(dir(fullfile(stage, '*.txt'))) == 48);
S = process_LocomToSleep(timeInterval=1, segmentMode=1, Teliminate=0, startZT=0, ...
    dayLength=24, deadThresh=12, minDays=3, bridgeWakeSeconds=1, ...
    folderPath=stage, selectAll=true);
fprintf('flies retained by the pipeline: %d of 48\n', numel(S.includedFiles));
types = ["sleep","wake"];
for ti = 1:2
    [B, meta] = extractBoutsForModelCompetition(S, 'entire', char(types(ti)), Verbose=false);
    O = meta.OnsetBins;
    for i = 1:size(O, 1)
        o = O(i, ~isnan(O(i, :)));
        assert(all(diff(o) > 0), 'bouts of fly %d not in time order', i);
    end
    nb = sum(~isnan(B), 2);
    keep = nb >= 40;
    names = string(S.includedFiles(:));
    fprintf('\n===== %s: %d flies in, %d with >= 40 bouts; bouts per kept fly %d-%d (median %g); min bout %g s\n', ...
        types(ti), numel(nb), sum(keep), min(nb(keep)), max(nb(keep)), median(nb(keep)), meta.MinObservedSeconds);
    if any(~keep), fprintf('  dropped: %s\n', strjoin(names(~keep) + " (" + nb(~keep) + ")", ', ')); end
    tic
    R = BoutIrreversibility(B(keep, :), Lags=[1 2 4 8], NumSurrogates=1000, Null="both", ...
        Seed=1, MinPairs=20, Plot=true, FlyNames=names(keep).');
    fprintf('  BoutIrreversibility: %.0f s\n', toc);
    disp(R.Group)
    exportgraphics(R.Figure, fullfile(out, "per0dd_" + types(ti) + "_group.png"), Resolution=150);
    close(R.Figure); R = rmfield(R, 'Figure');
    res.(types(ti)) = R; res.(types(ti) + "_bouts") = B(keep, :); res.(types(ti) + "_meta") = meta;
    writetable(R.PerFly, fullfile(out, "per0dd_" + types(ti) + "_perfly.csv"));
end
% primary family: lag 1, IAAFT, UpFraction and Bowley, sleep and wake; Holm
lab = ["sleep UpFraction","sleep Bowley","wake UpFraction","wake Bowley"];
p = [res.sleep.Group.UpFraction_P_IAAFT(1), res.sleep.Group.Bowley_P_IAAFT(1), ...
     res.wake.Group.UpFraction_P_IAAFT(1),  res.wake.Group.Bowley_P_IAAFT(1)];
[ps, ix] = sort(p); m = numel(p); adj = zeros(1, m); run = 0;
for k = 1:m, run = max(run, min(1, (m - k + 1) * ps(k))); adj(ix(k)) = run; end
fprintf('\n===== PRIMARY (lag 1, IAAFT, Holm across 4)\n');
for k = 1:m, fprintf('  %-18s p = %.4f   Holm-adjusted = %.4f\n', lab(k), p(k), adj(k)); end
res.primary = table(lab.', p.', adj.', 'VariableNames', {'Test','P','HolmP'});
res.Settings.stage = stage; res.Settings.when = string(datetime('now'));
save(fullfile(out, 'per0dd_results.mat'), 'res', '-v7.3');
