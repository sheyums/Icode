function G = gate_20r(stageDir, Teliminate, codeDir, minDays, tag)
% GATE_20R  Step 0 of 20r: re-extract per0 DD bouts and check them against the
% recorded statistics (PRESPEC "2026-10-31 20r gate amendment: stored bout files
% not found"). Run BY NAME from the Icode folder, with analysis_20r on the path.
%
%   stageDir   : folder holding ONLY the 25 per0 DD count files
%                (process_LocomToSleep with selectAll=true reads every .txt there).
%   Teliminate : hours trimmed from the start by process_LocomToSleep (0 for the
%                pre-cut cutdata_4.8.2025 files, per the 2026-10-31 PRESPEC entry).
%   codeDir    : folder holding the extraction pipeline (default: the current
%                Rockefeller folder; amendment 3 passes a git-archive of 04e0508^).
%   minDays    : process_LocomToSleep minDays (default 3).
%   tag        : suffix for the output files (default "T<Teliminate>").
%
% Writes analysis_20r/per0DD_bouts_20r_<tag>.mat and gate_20r_report_<tag>.txt
% (the first, failed run, with Teliminate 24, wrote the names without a suffix).

if nargin < 3 || isempty(codeDir)
    codeDir = 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts';
end
if nargin < 4 || isempty(minDays), minDays = 3; end
if nargin < 5 || isempty(tag), tag = sprintf('T%g', Teliminate); end

assert(isequal(fileparts(which('BoutHazard')), pwd), ...
    'Run from the Icode folder: BoutHazard resolves to %s', which('BoutHazard'));

rock = codeDir;
addpath(rock, '-begin');                 % ahead of any copy already on the saved path
cleanupPath = onCleanup(@() rmpath(rock));  % later steps use whatever the path was before
for f = {'process_LocomToSleep','sleeptrace','ironout','segment_sleep_phases', ...
         'extractBoutsForModelCompetition'}
    assert(isequal(fileparts(which(f{1})), rock), '%s resolves to %s', f{1}, which(f{1}));
end
% Every file the pipeline calls, from the folder actually used, with which -all.
deps = matlab.codetools.requiredFilesAndProducts( ...
    {which('process_LocomToSleep'), which('extractBoutsForModelCompetition')});
whichAll = {};
for k = 1:numel(deps)
    [~, nm] = fileparts(deps{k});
    if startsWith(deps{k}, rock)
        assert(isequal(fileparts(which(nm)), rock), '%s resolves to %s', nm, which(nm));
    end
    wa = which(nm, '-all'); if ischar(wa), wa = {wa}; end
    whichAll{end+1} = sprintf('%s:\n    %s', nm, strjoin(wa, sprintf('\n    '))); %#ok<AGROW>
end
for f = {'FitHyperexponentialMLE','FitTruncatedDiscreteMLE'}
    assert(isequal(fileparts(which(f{1})), pwd), '%s resolves to %s', f{1}, which(f{1}));
end

here = fileparts(mfilename('fullpath'));

% ------------------------------------------------------------- extraction
S = process_LocomToSleep(timeInterval=1, bridgeWakeSeconds=1, Teliminate=Teliminate, ...
    deadThresh=12, minDays=minDays, folderPath=stageDir, selectAll=true);
[Wm, metaW] = extractBoutsForModelCompetition(S, 'entire', 'wake');
[Sm, metaS] = extractBoutsForModelCompetition(S, 'entire', 'sleep');

dt    = S.timeInterval;
xminW = (floor(S.bridgeWakeSeconds/dt) + 1) * dt;
xminS = ceil(300/dt) * dt;

nFly   = numel(S.includedFiles);
w      = Wm(isfinite(Wm) & Wm >= xminW);
s      = Sm(isfinite(Sm) & Sm >= xminS);
wPerFly = sum(isfinite(Wm) & Wm >= xminW, 2);
sPerFly = sum(isfinite(Sm) & Sm >= xminS, 2);
belowW  = nnz(isfinite(Wm) & Wm < xminW);
belowS  = nnz(isfinite(Sm) & Sm < xminS);
flyDaysBouts = (sum(w) + sum(s)) / 86400;

% Edge diagnostic: is the record's first sample, and its last, inside a kept bout?
traceLen = zeros(nFly, 1); firstKept = false(nFly, 1); lastKept = false(nFly, 1);
for i = 1:nFly
    traceLen(i) = numel(S.(sprintf('fly%d', i)).entire);
    on  = [metaW.OnsetBins(i, :), metaS.OnsetBins(i, :)];
    dur = [Wm(i, :), Sm(i, :)] / dt;
    ok  = isfinite(on);
    firstKept(i) = any(on(ok) == 1);
    lastKept(i)  = any(on(ok) + dur(ok) - 1 == traceLen(i));
end

res = struct();
res.flies = nFly == 25;
res.wake  = numel(w) == 3989;
res.sleep = numel(s) == 3987;

% ------------------------------------------------------- 16-bin hazard table
fc = [1.41 2.45 3.87 6.71 12.0 21.2 36.7 63.1 109.3 190.1 330.8 574.5 997.7 1733 3012];
fd = [2 2 2 2 1 1 1 1 1 1 1 1 1 0 0];                  % decimals printed in FINDINGS
fr = [3989 3943 3909 3840 3755 3643 3545 3420 3260 3023 2623 2030 1186 455 104];
fe = [46 34 69 85 112 98 125 160 237 400 593 844 731 351 97];
H16 = BoutHazard(w, xminW, SamplingInterval=dt, NumBins=16, Spacing="log", ...
    MinAtRisk=10, Alpha=0.05);
nk = numel(H16.AtRisk);
if nk == numel(fr)
    res.centres = all(round(H16.BinCenters .* 10.^fd) == round(fc .* 10.^fd));
    res.atRisk  = isequal(H16.AtRisk, fr);
    res.events  = isequal(H16.Events, fe);
else
    res.centres = false; res.atRisk = false; res.events = false;
end

% Amendment 3: items (3)-(4) are computed only for a configuration whose counts
% and 16-bin table match exactly; otherwise the run is descriptive only.
matched = res.flies && res.wake && res.sleep && res.atRisk && res.events;
nb = [12 16 20 24]; rr = nan(1, 4); tr = nan(1, 4); pk = nan(1, 4);
tau = [NaN NaN]; q = [NaN NaN]; f2 = struct('LogLik', NaN);
if matched
    % --------------------------------------------- two-exponential null fit
    F  = FitHyperexponentialMLE(w, xminW, MaxComponents=2, DistributionType="discrete", ...
        SamplingInterval=dt, RandomSeed=1, Verbose=false);
    f2 = F.AllFits(2);
    tau = f2.Tau(:)'; q = f2.WeightsObserved(:)';
    res.tau  = isequal(round(tau, 1), [8.4 661.9]);
    res.q    = isequal(round(q, 3), [0.079 0.921]);
    res.logL = abs(f2.LogLik - (-29328.2460)) <= 1e-4;

    % --------------------------------------------------- RiseRatio, 4 counts
    for i = 1:4
        Hi = BoutHazard(w, xminW, SamplingInterval=dt, NumBins=nb(i), Spacing="log", ...
            MinAtRisk=10, Alpha=0.05);
        rr(i) = Hi.RiseRatio; tr(i) = Hi.BinCenters(Hi.TroughIndex); pk(i) = Hi.BinCenters(Hi.PeakIndex);
    end
    res.rise = isequal(round(rr, 3), [1.405 1.388 1.376 1.430]);
else
    res.tau = false; res.q = false; res.logL = false; res.rise = false;
end

% -------------------------------------------------------------- save bouts
files = S.includedFiles;
matPath = fullfile(here, sprintf('per0DD_bouts_20r_%s.mat', tag));
settings = struct('timeInterval', dt, 'bridgeWakeSeconds', S.bridgeWakeSeconds, ...
    'Teliminate', Teliminate, 'deadThresh', 12, 'minDays', minDays, 'phase', 'entire', ...
    'xminWake', xminW, 'xminSleep', xminS, 'codeDir', codeDir);
WakeBouts = Wm; SleepBouts = Sm; %#ok<NASGU>
save(matPath, 'WakeBouts', 'SleepBouts', 'files', 'settings', 'metaW', 'metaS', '-v7');

% ------------------------------------------------------------------- report
names = fieldnames(res); pass = all(struct2array(res));
L = {};
L{end+1} = sprintf('20r gate, %s, MATLAB %s', char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
L{end+1} = sprintf('stage dir: %s; Teliminate %g h; minDays %g; tag %s', stageDir, Teliminate, minDays, tag);
L{end+1} = sprintf('code dir: %s', codeDir);
L{end+1} = sprintf('pipeline files (%d), which -all, first entry is the one used:', numel(deps));
L = [L, whichAll];
L{end+1} = sprintf('flies %d; wake bouts >= %g s: %d (below xmin: %d); sleep bouts >= %g s: %d (below: %d)', ...
    nFly, xminW, numel(w), belowW, xminS, numel(s), belowS);
L{end+1} = sprintf('fly-days covered by bouts: %.2f', flyDaysBouts);
L{end+1} = sprintf(['per fly: wake %d-%d (recorded 81-248), sleep %d-%d (recorded 80-249); ' ...
    'first bout kept in %d flies, last bout kept in %d'], min(wPerFly), max(wPerFly), ...
    min(sPerFly), max(sPerFly), nnz(firstKept), nnz(lastKept));
L{end+1} = 'per fly (file, wake, sleep, record s, first kept, last kept):';
for i = 1:nFly
    L{end+1} = sprintf('  %2d %-70s %4d %4d %7d %d %d', i, files{i}, wPerFly(i), sPerFly(i), ...
        traceLen(i), firstKept(i), lastKept(i)); %#ok<AGROW>
end
L{end+1} = sprintf('16 bins: %d kept bins', nk);
L{end+1} = '  centre    atRisk  events   (FINDINGS: centre atRisk events)';
for j = 1:nk
    if j <= numel(fr), ref = sprintf('%g %d %d', fc(j), fr(j), fe(j)); else, ref = '-'; end
    L{end+1} = sprintf('  %8.2f  %6d  %6d   (%s)', H16.BinCenters(j), H16.AtRisk(j), H16.Events(j), ref); %#ok<AGROW>
end
L{end+1} = sprintf('null K=2: tau %.4f %.4f s, q %.5f %.5f, logL %.4f (recorded -29328.2460)', ...
    tau(1), tau(2), q(1), q(2), f2.LogLik);
for i = 1:4
    L{end+1} = sprintf('bins %2d: RiseRatio %.4f, trough %.1f s, peak %.1f s', nb(i), rr(i), tr(i), pk(i)); %#ok<AGROW>
end
for i = 1:numel(names)
    L{end+1} = sprintf('check %-8s %s', names{i}, string(res.(names{i}))); %#ok<AGROW>
end
L{end+1} = sprintf('GATE %s', ternary(pass, 'PASS', 'FAIL'));
L{end+1} = sprintf('saved %s', matPath);
txt = strjoin(L, newline);
fprintf('%s\n', txt);
fid = fopen(fullfile(here, sprintf('gate_20r_report_%s.txt', tag)), 'w'); fprintf(fid, '%s\n', txt); fclose(fid);

G = struct('res', res, 'pass', pass, 'nWake', numel(w), 'nSleep', numel(s), ...
    'H16', H16, 'NullFit', f2, 'RiseRatio', rr, 'files', {files});
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
