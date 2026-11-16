function [R, filePath] = LoadingCheck(varargin)
%LOADINGCHECK  Quick tube check on the first hours after loading a multibeam monitor.
%
%   R = LOADINGCHECK(FILEPATH) reads the first 2 h of a raw multibeam file (3 header lines; 64 positions and a
%   time column per line; the all -1 lines are skipped) and prints a RE-SEAT list: tubes whose beam readings
%   suggest that food, the cotton plug or debris is blocking a beam, that part of the tube is obstructed, or
%   that the fly hardly moves. Two hours come out of the first 24 h that are discarded anyway (entrainment).
%
%   Flags (per channel, first 'Hours' only):
%     BLOCKED BEAM k (end | mid-tube)  >= 10 jumps of >= 4 beams between consecutive readings, one value k
%                       receives >= 30% of them, and >= 20% of the visits to k that start with such a jump end
%                       by returning to the reading before it. A blocked beam's reading overrides the fly's
%                       position, so the reading jumps to k and back -- something a walking fly cannot do in
%                       0.4 s. An END beam (1, 2, 15, 16) points to food or the plug; mid-tube to debris,
%                       a droplet or a mark on the tube.
%     NARROW RANGE      the fly moves (>= 50 position changes) but spans < 8 beams: part of the tube obstructed.
%     FROZEN AT END BEAM k  one value in >= 99% of readings, at beam 1, 2, 15 or 16: the fly may be HIDDEN behind a
%                       blocked end beam (food or plug), or dead. In the 2026 pilot, 14 fly-weeks read a constant 16
%                       for the whole week while the fly was alive (active again later): re-seat these first.
%                       FROZEN AT BEAM k (mid-tube) and FROZEN AT 0 (no beam broken; dead flies often read 0) likewise.
%     LOW ACTIVITY      < 20 position changes: a dead or stuck fly, or a blocked tube.
%   A fly resting on one beam is NOT flagged by "reads one beam a lot"; only the jump signature counts.
%   These are the beam_hmm hub criteria (PRESPEC P1-C) with a minimum jump count for short recordings.
%
%   BLOCKED vs FROZEN, and WHAT TO DO (also printed after the RE-SEAT list, for the flags present):
%     BLOCKED    blocked SOMETIMES: the fly is visible most of the time, with impossible jumps to beam k and
%                back. Data usable after repair (RawToActivity1s / dead-hub rule). Re-seat if easy.
%     FROZEN     blocked ALL THE TIME (or a dead fly): one reading >= 99% of the time; a living fly can be
%                completely hidden. Data LOST while it lasts. Re-seat NOW: check that the plug or food is
%                not covering the end beam.
%     NARROW     part of the tube obstructed: check how far the food and plug reach into the tube.
%     LOW ACT    dead or stuck fly, or a blocked tube: check the fly; re-seat or replace.
%     Both BLOCKED and FROZEN (e.g. ch32 in the pilot): mostly hidden, occasional glimpses -- treat as FROZEN.
%
%   LIMITS (2026 pilot, 9 weekly files x 16 tubes): the 2-h screen flagged all 14 live flies hidden behind a
%   blocked beam 16, but only about half of the intermittent blockages found later in the full week (~70% with
%   any flag); some blockages develop during the week. Always run the full-file QC as well (RawToActivity1s
%   or BeamFingerprint) before scoring sleep.
%
%   R = LOADINGCHECK() or R = LOADINGCHECK('Channels', 17:32, ...) -- no file given -- opens a file dialog
%   (uigetfile; it starts in the folder of the last file picked). With 'Previous' and no file, the file of
%   the Previous check is used instead (see 'Previous').
%   [R, FILE] = LOADINGCHECK(...) also returns the full path of the file read.
%
%   OUTPUT: a line per channel in the Command Window (OK or its flags, with the numbers behind them), the
%   RE-SEAT list, a table R, and figures of 16 channels each (4 x 4 panels; 64 channels -> 4 figures, which
%   MATLAB R2025a+ shows as tabs of one window): the position trace over the window (0 = beyond beam 1 ... 16),
%   the suspected blocked beam as a red line, the title red when flagged.
%
%   R (table; a struct of columns in Octave): Channel, Flags (full wording), Changes (position changes),
%   BigJumps (>= 4 beams), HubValue/HubShare/HubReturns (the most frequent jump target and the two hub
%   criteria), BeamsSpanned. R.Properties.UserData holds File (full path), StartHour and Hours (Octave: the
%   fields R.File, R.StartHour, R.Hours); R.Properties.Description is the full path as well.
%
%   NAME-VALUE:
%     'Hours'      length of the window checked (default 2 -- RECOMMENDED). 1 h is acceptable for a first quick
%                  look: in the pilot (144 fly-weeks), the first 1 h and the first 2 h flagged the same blocked
%                  beams (18 + 1) and the same 14 hidden flies; 1 h added one NARROW flag. But the thresholds
%                  are fixed counts, so in a quieter hour a fly gives fewer jumps (BLOCKED missed) and a resting
%                  fly fills more of the window (FROZEN / LOW ACTIVITY on a sleeping fly). A warning is printed
%                  for windows shorter than 2 h.
%                  Check soon after loading or re-seating, while handling keeps the flies active; a window in
%                  the midday siesta or at night flags sleeping flies as LOW ACTIVITY or FROZEN.
%     'StartHour'  where the window starts, in hours after the file's first reading (default 0). Use it to
%                  check tubes AFTER a re-seat in the same, still-running file: re-seated at ~2.3 h into the
%                  file -> 'StartHour', 2.5 (leave a few minutes for handling).
%     'Previous'   the R returned by an earlier LoadingCheck (e.g. before the re-seat). Each channel is then
%                  labelled RESOLVED (flagged before, OK now), STILL FLAGGED, NEW FLAG (OK before, flagged
%                  now -- re-seating can disturb a neighbour) or OK (OK before and now); the printout lists
%                  before -> after, R gains PreviousFlags and Change, and the panel titles show the change.
%                  With no file given, the Previous check's file is read again (R2 = LoadingCheck('StartHour',
%                  2.5, 'Previous', R1)). If a file IS given and differs from the Previous one, a warning names
%                  both (fine after starting a new file; otherwise a wrong file). A Previous from a version
%                  without the stored path is compared as before, unchecked.
%     'Channels' (default 1:64), 'Plot' (default true), 'SaveFile' (PNG path; default none; with more than
%     16 channels one PNG per figure, named <SaveFile>_ch1-16.png, ...; with a SaveFile the figures get a
%     fixed size, so the PNGs are reproducible, and do not dock as tabs),
%     'Visible' ('on'/'off'). With 'Visible','off' the figures are closed at the end (after saving).
%     For headless runs end the batch command with: close all; exit
%
%   AFTER RE-SEATING: did it work?
%     R1 = LoadingCheck(f, 'Channels', 17:32);                       % first 2 h: RE-SEAT list
%     ... re-seat the flagged tubes at, say, 2.3 h into the recording; keep recording >= 2 h more ...
%     R2 = LoadingCheck(f, 'Channels', 17:32, 'StartHour', 2.5, 'Previous', R1);
%     R2 = LoadingCheck('Channels', 17:32, 'StartHour', 2.5, 'Previous', R1);   % same, f taken from R1
%     RESOLVED only means "not seen in this window": about half of the intermittent blockages in the pilot
%     appeared later in the week, so the full-file QC is still needed. (Or start a new file after re-seating
%     and check its first 2 h with 'Previous', R1.)
%
%   EXAMPLE
%     R = LoadingCheck('raw_positions_all64_week2.txt', 'Channels', 17:32, 'SaveFile', 'week2_loading.png');
%
%   TYPICAL SEQUENCE (functions in beam_hmm/matlab: cd there or addpath it; one weekly raw file at a time)
%     1  After loading the tubes, record ~2 h, then screen the tubes and re-seat the flagged ones:
%          R = LoadingCheck('raw_positions_all64_week2.txt', 'Channels', 17:32);
%        After re-seating (here at ~2.3 h), record >= 2 h more and confirm:
%          R2 = LoadingCheck('raw_positions_all64_week2.txt', 'Channels', 17:32, 'StartHour', 2.5, 'Previous', R);
%     2  After the recording, look at the whole file (QC table, channel status, four figures):
%          F = BeamFingerprint('raw_positions_all64_week2.txt', 'Channels', 17:32, 'SaveFolder', 'fp_week2');
%     3  Convert to activity per second (repairs + QC), one .txt per channel:
%          A = RawToActivity1s('raw_positions_all64_week2.txt', 'Channels', 17:32, 'OutputFolder', 'week2_1s');
%     4  Score sleep with the lab pipeline, unchanged (one experiment per folder; set startZT -- nothing checks it):
%          S = process_LocomToSleep(timeInterval=1, folderPath='week2_1s', selectAll=true, startZT=0);
%     Optional: EstimatePhaseStamp to check the phase (startZT) from behaviour; chi2p for the period in DD.
%     Once per installation: test_BeamFingerprint and test_RawToActivity1s on the week-1 raw file (golden tests).
%   Python twin: beam_hmm/labconv.py loading_check (identical flags; validated against the full-week QC).
%   See also RAWTOACTIVITY1S, BEAMFINGERPRINT, PROCESS_LOCOMTOSLEEP.

persistent lastDir                                               % folder of the last file picked in the dialog
pNames = {'Hours', 'Channels', 'Plot', 'SaveFile', 'Visible', 'StartHour', 'Previous'};
isTxt = nargin >= 1 && (ischar(varargin{1}) || isstring(varargin{1}));
if isTxt, a1 = char(varargin{1}); else, a1 = ''; end
isName = ~isempty(a1) && any(strncmpi(a1, pNames, numel(a1))) && ~exist(a1, 'file');  % a parameter name,
if isTxt && ~isempty(a1) && ~isName                                 % or a prefix of one ('Chan'), as inputParser
    filePath = a1; args = varargin(2:end);                       % LoadingCheck(f, ...)   accepts those too
elseif nargin >= 1 && (isempty(varargin{1}) || (isTxt && isempty(a1)))   % [], '' or ""
    filePath = ''; args = varargin(2:end);                       % LoadingCheck([], ...) or LoadingCheck('', ...)
else
    filePath = ''; args = varargin;                              % LoadingCheck() or LoadingCheck('Channels', ...)
end
p = inputParser;
p.addParameter('Hours', 2); p.addParameter('Channels', 1:64);
p.addParameter('Plot', true); p.addParameter('SaveFile', ''); p.addParameter('Visible', 'on');
p.addParameter('StartHour', 0, @(v) isnumeric(v) && isscalar(v) && v >= 0); p.addParameter('Previous', []);
p.parse(args{:}); o = p.Results;
prevFile = previousFile(o.Previous);                             % '' if Previous is absent or predates this field
if isempty(filePath) && ~isempty(prevFile)
    filePath = prevFile;                                         % re-check the file the Previous check read
    fprintf('LoadingCheck: no file given; using the Previous check''s file\n');
elseif isempty(filePath)
    if isempty(lastDir) || ~exist(lastDir, 'dir'), lastDir = pwd; end
    [fn, fd] = uigetfile({'*.txt', 'Raw multibeam files (*.txt)'; '*.*', 'All files (*.*)'}, ...
                         'LoadingCheck: pick a raw multibeam file', [lastDir filesep]);
    if isequal(fn, 0), error('LoadingCheck:NoFileChosen', 'No file chosen (dialog cancelled).'); end
    filePath = fullfile(fd, fn); lastDir = fd;
end
fid = fopen(filePath, 'r');
if fid < 0, error('LoadingCheck:NoFile', 'Cannot open %s', filePath); end
filePath = fullPathOf(fopen(fid));                               % absolute, '.' and '..' resolved
if ~isempty(prevFile) && ~sameFile(filePath, prevFile)
    warning('LoadingCheck:DifferentFile', ['Previous came from a DIFFERENT file:\n  previous: %s\n  now:      %s\n' ...
            'Fine if you started a new file after re-seating; otherwise check the file name.'], prevFile, filePath);
end
for i = 1:3, fgetl(fid); end
nLines = 2 * ceil((o.StartHour + o.Hours) * 3600 / 0.2) + 2000;        % ample for the window at <= 5 Hz, -1 lines and CRLF line ends included (Octave counts CRLF twice)
C = textscan(fid, [repmat('%d8', 1, 64) '%f'], nLines, 'Delimiter', '\t', 'EmptyValue', NaN, 'CollectOutput', true);
fclose(fid);
A = C{1}; t = C{2}; keep = ~isnan(t); A = A(keep, :); t = t(keep);
rel = t - t(1); m = rel >= o.StartHour * 3600 & rel < (o.StartHour + o.Hours) * 3600;
if ~any(m)
    error('LoadingCheck:NoData', 'No readings %.2f-%.2f h after the file start (the file holds %.2f h).', ...
          o.StartHour, o.StartHour + o.Hours, rel(end) / 3600);
end
tw = rel(m); A = double(A(m, :));
if (tw(end) - o.StartHour * 3600) < 0.95 * o.Hours * 3600
    warning('LoadingCheck:Short', 'Only %.2f h of readings in the window; checked what there is.', ...
            (tw(end) - o.StartHour * 3600) / 3600);
end
if o.Hours < 2
    warning('LoadingCheck:ShortWindow', ['%.1f-h window: thresholds are fixed counts, so BLOCKED can be missed and ' ...
            'FROZEN / LOW ACTIVITY can fire on a resting fly more easily than in the recommended 2 h.'], o.Hours);
end
ch = o.Channels(:)'; nC = numel(ch);
[nchg, nbig, hubV, hubSh, hubRet, span] = deal(nan(nC, 1)); flags = repmat({''}, nC, 1);
for k = 1:nC
    x = A(:, ch(k)); d = diff(x); nchg(k) = sum(d ~= 0); big = find(abs(d) >= 4) + 1; nbig(k) = numel(big);
    span(k) = max(x) - min(x); f = {};
    if ~isempty(big)
        [v, ~, j] = unique(x(big)); cnt = accumarray(j, 1); [mx, im] = max(cnt); hubV(k) = v(im); hubSh(k) = mx / numel(big);
        st = big(x(big) == hubV(k)); r = false(numel(st), 1);
        for s = 1:numel(st)
            i = st(s); q = i;
            while q + 1 <= numel(x) && x(q + 1) == hubV(k), q = q + 1; end
            r(s) = (q + 1 <= numel(x)) && x(q + 1) == x(i - 1);
        end
        hubRet(k) = mean(r);
        if nbig(k) >= 10 && hubSh(k) >= 0.30 && hubRet(k) >= 0.20
            if any(hubV(k) == [1 2 15 16]), w = 'end'; else, w = 'mid-tube'; end
            f{end + 1} = sprintf('BLOCKED BEAM %d (%s)', hubV(k), w);
        end
    end
    if nchg(k) >= 50 && span(k) < 8, f{end + 1} = 'NARROW RANGE'; end
    [vals, ~, jj] = unique(x); cnts = accumarray(jj, 1); [mxc, imc] = max(cnts); fv = vals(imc);
    if mxc / numel(x) >= 0.99
        if fv == 0, f{end + 1} = 'FROZEN AT 0 (no beam broken)';
        elseif any(fv == [1 2 15 16]), f{end + 1} = sprintf('FROZEN AT END BEAM %d', fv);
        else, f{end + 1} = sprintf('FROZEN AT BEAM %d', fv); end
    end
    if nchg(k) < 20, f{end + 1} = 'LOW ACTIVITY'; end
    flags{k} = strjoin(f, '; ');
end
bad = find(~cellfun(@isempty, flags));
if o.StartHour == 0
    fprintf('LoadingCheck: %s, first %.1f h, %d channels\n', filePath, o.Hours, nC);
else
    fprintf('LoadingCheck: %s, %.2f-%.2f h after the file start, %d channels\n', filePath, o.StartHour, ...
            o.StartHour + o.Hours, nC);
end
% ---- comparison with an earlier check ('Previous') ----
hasPrev = ~isempty(o.Previous); prevF = repmat({''}, nC, 1); change = repmat({''}, nC, 1);
if hasPrev
    pc = double(o.Previous.Channel(:)); pf = o.Previous.Flags; if ischar(pf), pf = {pf}; end
    for k = 1:nC
        i = find(pc == ch(k), 1);
        if isempty(i), change{k} = 'NOT IN PREVIOUS'; prevF{k} = ''; continue; end
        prevF{k} = pf{i}; b = ~isempty(prevF{k}); a = ~isempty(flags{k});
        if b && ~a, change{k} = 'RESOLVED';
        elseif b && a, change{k} = 'STILL FLAGGED';
        elseif ~b && a, change{k} = 'NEW FLAG';
        else, change{k} = 'OK'; end
    end
end
fprintf('  %-6s %8s %9s %5s %6s %8s %6s   %s\n', 'chan', 'changes', 'bigjumps', 'hub', 'share', 'returns', 'span', 'verdict');
for k = 1:nC
    v = flags{k}; if isempty(v), v = 'OK'; end
    fprintf('  ch%-4d %8d %9d %5s %6.2f %8.2f %6d   %s\n', ch(k), nchg(k), nbig(k), num2str(hubV(k)), ...
            hubSh(k), hubRet(k), span(k), v);
end
if isempty(bad)
    fprintf('  no tube flagged\n');
else
    fprintf('  RE-SEAT / CHECK %d tube(s):\n', numel(bad));
    for k = bad(:)', fprintf('    ch%-3d %s\n', ch(k), flags{k}); end
    allf = strjoin(flags(bad)', ' ');
    fprintf('  WHAT TO DO:\n');
    if ~isempty(strfind(allf, 'FROZEN'))
        fprintf('    FROZEN   blocked all the time (or dead fly): a live fly can be HIDDEN, data lost while it lasts.\n');
        fprintf('             Re-seat NOW; check that the plug or food is not covering the end beam.\n');
    end
    if ~isempty(strfind(allf, 'BLOCKED'))
        fprintf('    BLOCKED  blocked sometimes: fly visible most of the time; data usable after repair. Re-seat if easy.\n');
    end
    if ~isempty(strfind(allf, 'NARROW')), fprintf('    NARROW   part of the tube obstructed: check how far food and plug reach into the tube.\n'); end
    if ~isempty(strfind(allf, 'LOW ACTIVITY')), fprintf('    LOW ACT  dead or stuck fly, or a blocked tube: check the fly; re-seat or replace.\n'); end
    fprintf('  (A 2-h screen misses blockages that develop later: run the full-file QC before scoring sleep.)\n');
end
if hasPrev
    fprintf('  COMPARED WITH THE PREVIOUS CHECK (before -> now):\n');
    for k = 1:nC
        if strcmp(change{k}, 'OK'), continue; end
        b = prevF{k}; if isempty(b), b = 'OK'; end
        if strcmp(change{k}, 'NOT IN PREVIOUS'), b = '(not checked)'; end
        a = flags{k}; if isempty(a), a = 'OK'; end
        fprintf('    ch%-3d %-15s %s  ->  %s\n', ch(k), change{k}, b, a);
    end
    nRes = sum(strcmp(change, 'RESOLVED')); nSt = sum(strcmp(change, 'STILL FLAGGED'));
    nNew = sum(strcmp(change, 'NEW FLAG')); nOK = sum(strcmp(change, 'OK'));
    nNot = sum(strcmp(change, 'NOT IN PREVIOUS'));
    fprintf('    %d resolved, %d still flagged, %d new, %d OK before and now', nRes, nSt, nNew, nOK);
    if nNot > 0, fprintf(', %d not in the previous check', nNot); end
    fprintf('\n');
    if nSt > 0, fprintf('    STILL FLAGGED: re-seat again or look at the tube and fly (food/plug position, fly stuck or dead).\n'); end
    if nNew > 0, fprintf('    NEW FLAG: check whether handling disturbed this tube or its fly.\n'); end
    if nRes > 0, fprintf('    RESOLVED means not seen in this window only; run the full-file QC after the recording.\n'); end
end
if o.Plot
    tt = (tw - o.StartHour * 3600) / 60;                           % minutes since the window start
    perFig = 16; nPage = ceil(nC / perFig);                        % 16 channels per figure (4 x 4), one figure
    [~, fstem] = fileparts(filePath);                              % per 16 -- R2025a+ opens them as tabs
    for pg = 1:nPage
    kk = (pg - 1) * perFig + 1:min(pg * perFig, nC); nOn = numel(kk);
    nCol = min(4, nOn); nRow = ceil(nOn / nCol);
    figArgs = {'Color', 'w', 'Visible', o.Visible, 'Name', sprintf('LoadingCheck ch%d-%d', ch(kk(1)), ch(kk(end))), ...
               'NumberTitle', 'off'};
    if ~isempty(o.SaveFile) || strcmpi(o.Visible, 'off') || isOctave   % saved: a fixed size, so the PNG is
        figArgs = [figArgs {'Position', [50 50 260 * nCol 210 * nRow + 70]}];   % reproducible; shown, unsaved
    end                                                            % MATLAB figures keep the default size and dock as tabs
    fh = figure(figArgs{:});
    try, theme(fh, 'light'); catch, end                          % pin the theme before any axes (CLAUDE.md)
    try, sgtitle(fh, sprintf('%s   ch%d-%d', fstem, ch(kk(1)), ch(kk(end))), 'FontSize', 9, 'Interpreter', 'none'); catch, end
    for k = kk
        ax = subplot(nRow, nCol, k - kk(1) + 1, 'Parent', fh); hold(ax, 'on');
        if ~isempty(flags{k}) && ~isempty(strfind(flags{k}, 'BLOCKED'))
            plot(ax, [0 o.Hours * 60], [hubV(k) hubV(k)], '-', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5);
        end
        stairs(ax, tt, A(:, ch(k)), 'k', 'LineWidth', 0.5);
        ylim(ax, [-0.5 16.5]); xlim(ax, [0 o.Hours * 60]); set(ax, 'YTick', [0 8 16], 'FontSize', 7);
        if isempty(flags{k}), l2 = 'OK'; col = [0 0 0]; else, l2 = shortFlag(flags{k}); col = [0.85 0.1 0.1]; end
        if hasPrev
            switch change{k}
                case 'RESOLVED', l2 = 'RESOLVED'; col = [0 0.5 0];
                case 'STILL FLAGGED', l2 = ['STILL ' l2];
                case 'NEW FLAG', l2 = ['NEW ' l2];
                case 'NOT IN PREVIOUS', l2 = ['(new ch) ' l2];
            end
        end
        title(ax, {sprintf('ch%d', ch(k)), l2}, 'FontSize', 8, 'Color', col, 'Interpreter', 'none');
        j = k - kk(1) + 1;                                         % panel index within this figure
        if mod(j - 1, nCol) == 0, ylabel(ax, 'beam'); end
        if j > nOn - nCol
            if o.StartHour == 0, xlabel(ax, 'min'); else, xlabel(ax, sprintf('min after %.2f h', o.StartHour)); end
        end
    end
    if ~isempty(o.SaveFile)
        sf = char(o.SaveFile);
        if nPage > 1                                               % one PNG per figure: name_ch17-32.png, ...
            [sd, sn, se] = fileparts(sf); if isempty(se), se = '.png'; end
            sf = fullfile(sd, sprintf('%s_ch%d-%d%s', sn, ch(kk(1)), ch(kk(end)), se));
        end
        print(fh, sf, '-dpng', '-r150');
    end
    if strcmpi(o.Visible, 'off'), close(fh); end                   % never shown (saved or not): close it, so a
                                                                   % headless (-batch) MATLAB can exit and loops
                                                                   % do not pile up invisible figures
    end
end
names = {'Channel', 'Flags', 'Changes', 'BigJumps', 'HubValue', 'HubShare', 'HubReturns', 'BeamsSpanned'};
cols = {ch', flags, nchg, nbig, hubV, hubSh, hubRet, span};
if hasPrev, names = [names {'PreviousFlags', 'Change'}]; cols = [cols {prevF, change}]; end
info = struct('File', filePath, 'StartHour', o.StartHour, 'Hours', o.Hours);
if exist('table', 'file') || exist('table', 'builtin')
    R = table(cols{:}, 'VariableNames', names);
    try                                                            % R.Properties.UserData.File: full path read
        R.Properties.UserData = info; R.Properties.Description = filePath;
    catch                                                          % Octave table packages may lack these
        R = cell2struct(cols, names, 2); R.File = info.File; R.StartHour = info.StartHour; R.Hours = info.Hours;
    end
else
    R = cell2struct(cols, names, 2); R.File = info.File; R.StartHour = info.StartHour; R.Hours = info.Hours;
end
end

function f = previousFile(P)
% Full path stored by an earlier LoadingCheck, or '' (no Previous, or one made before the field existed).
f = '';
if isempty(P), return; end
if isa(P, 'table')                                                 % isa, not istable: Octave has no istable
    u = P.Properties.UserData; if isstruct(u) && isfield(u, 'File'), f = char(u.File); end
elseif isstruct(P) && isfield(P, 'File')
    f = char(P.File);
end
end

function f = fullPathOf(nm)
% Absolute path of an opened file. dir() gives the resolved folder ('.', '..' removed); a file that fopen
% found on the MATLAB path, not in pwd, is located with which(); the last resort prepends pwd.
f = '';
d = dir(nm);
if numel(d) == 1 && isfield(d, 'folder'), f = fullfile(d.folder, d.name); return; end
w = which(nm);
if ~isempty(w), f = w; return; end
f = nm;
if isempty(regexp(f, '^([A-Za-z]:[\\/]|[\\/~])', 'once')), f = fullfile(pwd, f); end
end

function tf = sameFile(a, b)
% Same file? Case-insensitive on Windows; slashes normalised.
a = strrep(a, '\', '/'); b = strrep(b, '\', '/');
if ispc, tf = strcmpi(a, b); else, tf = strcmp(a, b); end
end

function tf = isOctave
tf = exist('OCTAVE_VERSION', 'builtin') ~= 0;
end

function s = shortFlag(f)
% Compact panel label; the full wording is in the printout and in R.Flags.
s = regexprep(f, 'BLOCKED BEAM (\d+) \((end|mid-tube)\)', 'BLOCKED $1');
s = regexprep(s, 'FROZEN AT END BEAM (\d+)', 'FROZEN $1');
s = regexprep(s, 'FROZEN AT BEAM (\d+)', 'FROZEN $1');
s = strrep(s, 'FROZEN AT 0 (no beam broken)', 'FROZEN 0');
s = strrep(s, 'NARROW RANGE', 'NARROW'); s = strrep(s, 'LOW ACTIVITY', 'LOW ACT');
s = strrep(s, '; ', ' · ');
end
