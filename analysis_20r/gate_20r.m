function G = gate_20r(stageDir, Teliminate)
% GATE_20R  Step 0 of 20r: re-extract per0 DD bouts and check them against the
% recorded statistics (PRESPEC "2026-10-31 20r gate amendment: stored bout files
% not found"). Run BY NAME from the Icode folder, with analysis_20r on the path.
%
%   stageDir   : folder holding ONLY the 25 per0 DD count files
%                (process_LocomToSleep with selectAll=true reads every .txt there).
%   Teliminate : hours trimmed from the start by process_LocomToSleep (0 for the
%                pre-cut cutdata_4.8.2025 files, per the 2026-10-31 PRESPEC entry).
%
% Writes analysis_20r/per0DD_bouts_20r_T<Teliminate>.mat and gate_20r_report_T<Teliminate>.txt
% (the first, failed run, with Teliminate 24, wrote the names without the suffix).

assert(isequal(fileparts(which('BoutHazard')), pwd), ...
    'Run from the Icode folder: BoutHazard resolves to %s', which('BoutHazard'));

rock = 'C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts';
addpath(rock, '-end');
for f = {'process_LocomToSleep','sleeptrace','ironout','segment_sleep_phases', ...
         'extractBoutsForModelCompetition'}
    assert(isequal(fileparts(which(f{1})), rock), '%s resolves to %s', f{1}, which(f{1}));
end
for f = {'FitHyperexponentialMLE','FitTruncatedDiscreteMLE'}
    assert(isequal(fileparts(which(f{1})), pwd), '%s resolves to %s', f{1}, which(f{1}));
end

here = fileparts(mfilename('fullpath'));

% ------------------------------------------------------------- extraction
S = process_LocomToSleep(timeInterval=1, bridgeWakeSeconds=1, Teliminate=Teliminate, ...
    deadThresh=12, minDays=3, folderPath=stageDir, selectAll=true);
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

% ------------------------------------------------- two-exponential null fit
F  = FitHyperexponentialMLE(w, xminW, MaxComponents=2, DistributionType="discrete", ...
    SamplingInterval=dt, RandomSeed=1, Verbose=false);
f2 = F.AllFits(2);
tau = f2.Tau(:)'; q = f2.WeightsObserved(:)';
res.tau  = isequal(round(tau, 1), [8.4 661.9]);
res.q    = isequal(round(q, 3), [0.079 0.921]);
res.logL = abs(f2.LogLik - (-29328.2460)) <= 1e-4;

% ------------------------------------------------------- RiseRatio, 4 counts
nb = [12 16 20 24]; rr = nan(1, 4); tr = nan(1, 4); pk = nan(1, 4);
for i = 1:4
    Hi = BoutHazard(w, xminW, SamplingInterval=dt, NumBins=nb(i), Spacing="log", ...
        MinAtRisk=10, Alpha=0.05);
    rr(i) = Hi.RiseRatio; tr(i) = Hi.BinCenters(Hi.TroughIndex); pk(i) = Hi.BinCenters(Hi.PeakIndex);
end
res.rise = isequal(round(rr, 3), [1.405 1.388 1.376 1.430]);

% -------------------------------------------------------------- save bouts
files = S.includedFiles;
matPath = fullfile(here, sprintf('per0DD_bouts_20r_T%g.mat', Teliminate));
settings = struct('timeInterval', dt, 'bridgeWakeSeconds', S.bridgeWakeSeconds, ...
    'Teliminate', Teliminate, 'deadThresh', 12, 'minDays', 3, 'phase', 'entire', ...
    'xminWake', xminW, 'xminSleep', xminS);
WakeBouts = Wm; SleepBouts = Sm; %#ok<NASGU>
save(matPath, 'WakeBouts', 'SleepBouts', 'files', 'settings', 'metaW', 'metaS', '-v7');

% ------------------------------------------------------------------- report
names = fieldnames(res); pass = all(struct2array(res));
L = {};
L{end+1} = sprintf('20r gate, %s, MATLAB %s', char(datetime('now','Format','yyyy-MM-dd HH:mm')), version);
L{end+1} = sprintf('stage dir: %s; Teliminate %g h', stageDir, Teliminate);
L{end+1} = sprintf('flies %d; wake bouts >= %g s: %d (below xmin: %d); sleep bouts >= %g s: %d (below: %d)', ...
    nFly, xminW, numel(w), belowW, xminS, numel(s), belowS);
L{end+1} = sprintf('fly-days covered by bouts: %.2f', flyDaysBouts);
L{end+1} = 'per fly (file, wake, sleep):';
for i = 1:nFly
    L{end+1} = sprintf('  %2d %-50s %4d %4d', i, files{i}, wPerFly(i), sPerFly(i)); %#ok<AGROW>
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
fid = fopen(fullfile(here, sprintf('gate_20r_report_T%g.txt', Teliminate)), 'w'); fprintf(fid, '%s\n', txt); fclose(fid);

G = struct('res', res, 'pass', pass, 'nWake', numel(w), 'nSleep', numel(s), ...
    'H16', H16, 'NullFit', f2, 'RiseRatio', rr, 'files', {files});
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
