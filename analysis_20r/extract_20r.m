function B = extract_20r(stageDir, outMat)
% EXTRACT_20R  The 20r extractor (PRESPEC "20r amendment: stored per0 DD bouts
% lost"): the current Rockefeller pipeline, timeInterval 1, bridgeWakeSeconds 1,
% Teliminate 0, deadThresh 12, minDays 3, phase 'entire'. Bouts touching the
% record's first or last sample are dropped as censored (segment_sleep_phases),
% the convention Sheyum designated on 2026-10-31.
%
%   stageDir : folder holding ONLY the count files of one genotype
%   outMat   : .mat to write (wake and sleep bouts per fly, NaN-padded, seconds)
%
% Returns B with fields Wake, Sleep (fly x bout), xminWake, xminSleep, files,
% traceLen (samples per retained trace), metaW, metaS.

codeDir = 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts';
addpath(codeDir, '-begin');
cleanupPath = onCleanup(@() rmpath(codeDir));
for f = {'process_LocomToSleep','sleeptrace','ironout','segment_sleep_phases', ...
         'extractBoutsForModelCompetition'}
    assert(isequal(fileparts(which(f{1})), codeDir), '%s resolves to %s', f{1}, which(f{1}));
end

S = process_LocomToSleep(timeInterval=1, bridgeWakeSeconds=1, Teliminate=0, ...
    deadThresh=12, minDays=3, folderPath=stageDir, selectAll=true);
[Wm, metaW] = extractBoutsForModelCompetition(S, 'entire', 'wake');
[Sm, metaS] = extractBoutsForModelCompetition(S, 'entire', 'sleep');

dt = S.timeInterval;
B = struct();
B.Wake = Wm; B.Sleep = Sm;
B.xminWake  = (floor(S.bridgeWakeSeconds/dt) + 1) * dt;
B.xminSleep = ceil(300/dt) * dt;
B.files = S.includedFiles;
B.traceLen = cellfun(@(k) numel(S.(sprintf('fly%d', k)).entire), num2cell(1:numel(B.files)))';
B.metaW = metaW; B.metaS = metaS;
B.settings = struct('timeInterval', dt, 'bridgeWakeSeconds', S.bridgeWakeSeconds, ...
    'Teliminate', 0, 'deadThresh', 12, 'minDays', 3, 'phase', 'entire', 'codeDir', codeDir);
save(outMat, '-struct', 'B', '-v7');
end
