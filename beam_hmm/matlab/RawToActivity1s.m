function out = RawToActivity1s(filePath, varargin)
%RAWTOACTIVITY1S  Raw multibeam position file -> activity per second, in one step, with QC.
%
%   OUT = RAWTOACTIVITY1S(FILEPATH) reads one raw multibeam file (3 header lines; per line 64 tab-separated
%   positions and a time column; every other line is 64 x -1 with no time) and returns, per channel, the
%   number of BEAMS TRAVELLED in each second -- the same unit as the lab's multibeam_locomotion output, so
%   sleeptrace / process_LocomToSleep take it unchanged (activity <= 0 means no movement in that second).
%
%   It replaces PreProcess.m + multibeam_locomotion.m and offers two modes:
%
%   'beam_hmm' (default, recommended)
%       - the all -1 lines are DROPPED (the true sampling step, 0.4 s here, is measured and reported);
%       - spike filter: a single reading that jumps >= 4 beams and comes straight back (A-B-A) is set to A;
%       - dead-hub repair: if one value receives >= 30% of all >= 4-beam jumps AND >= 20% of the runs at it
%         that start with such a jump return to the value before it, that beam is treated as DEAD (as the
%         hardware already treats beam 15): every reading of it is replaced by the last other reading. This is
%         an intermittently blocked beam (debris, a droplet, food) whose reading overrides the fly's position;
%       - position 0 ("no beam broken": in these data always entered from beam 1, i.e. the fly beyond beam 1;
%         in a correctly oriented tube beam 1 is nearest the food, so 0 = fly at the food, very likely eating)
%         is KEPT, so crossing beam 1 counts as movement ('KeepBeyondBeam1', true). Set false to fold 0 into
%         the previous reading, as PreProcess does;
%       - second s = all position changes whose sample time t falls in [s, s+1), t from the file's time column.
%       Validated on synthetic flies with a planted blocked beam: sleep/quiet/active minutes agree 97-99% with
%       the artifact-free truth, and no artifact-free fly is changed (beam_hmm PRESPEC P1-C, R1 check).
%
%   'legacy'
%       Exactly PreProcess.m + multibeam_locomotion.m (binner = 5): -1, 0, and 1-not-preceded-by-2 replaced by
%       the previous value on ALL lines, |displacement| per line, sums over consecutive blocks of 5 lines.
%       Use it to reproduce older results. Blocked-beam flicker then counts as movement: in the 2026 pilot this
%       removed up to 261 min/day of sleep in affected fly-weeks.
%
%   QC (out.QC, and printed): per channel the status (EMPTY/DEAD, INCLUDED), hub beam and whether it was
%   declared dead, END-BEAM HUB (hub at beam 1, 2, 15 or 16: e.g. food blocking an end beam; the fly's
%   position at that end cannot be known), NARROW RANGE (< 8 beams spanned), and the share of time at 0.
%   QC is computed in both modes; repairs are applied only in 'beam_hmm' mode with 'Repair' true.
%
%   NAME-VALUE OPTIONS
%     'Mode'             'beam_hmm' (default) | 'legacy'
%     'Repair'           true (default): spike filter + dead-hub repair ('beam_hmm' mode only)
%     'KeepBeyondBeam1'  true (default): 0 is a position; false folds 0 into the previous reading
%     'Channels'         default 1:64
%     'OutputFolder'     if given, one file per channel is written: <FilePrefix><channel>.txt, one value per
%                        line ('%d \r\n', the lab's format). Existing files are NOT appended to:
%                        the call errors unless 'Overwrite' is true (multibeam_locomotion opened with 'a+',
%                        which silently doubled a trace when run twice).
%     'FilePrefix'       default: the raw file's name followed by '_ch'
%     'Overwrite'        default false
%     'AtZeroFolder'     if given, one file per channel with out.AtZero (see OUTPUT): <FilePrefix><channel>_at0.txt,
%                        one value per line ('%.4f'; NaN for a second with no reading). Must differ from
%                        OutputFolder: process_LocomToSleep(selectAll=true) reads every file in its folder.
%
%   OUTPUT
%     out.Counts    seconds x channels, beams travelled per second (double)
%     out.AtZero    seconds x channels, share of each second's readings at position 0 (0 .. 1; NaN if the second
%                   has no reading). In a correctly oriented tube (food next to beam 1) 0 means the fly is at the
%                   food: a FEEDING PROXY ("at the food", not ingestion). In a reversed tube it is the plug end.
%                   Computed after the same repairs as Counts, in every mode and whatever 'KeepBeyondBeam1' is.
%                   It does not change Counts or the sleep scoring: a fly still at 0 for >= 5 min is scored
%                   asleep; use AtZero to flag such bouts (pilot: <= 7.5% of any fly's sleep, median ~0%).
%     out.Channels, out.Mode, out.Settings, out.Info (header, step, irregular steps), out.QC
%     QC columns include ShareAtZero (share of readings at 0, raw) and EntriesToZeroPerDay (entries into 0 per
%     day of recording, after the repairs; pilot median ~155 per day in readable fly-weeks).
%
%   EXAMPLE
%     A = RawToActivity1s('raw_positions_all64_week2.txt', 'Channels', 17:32, 'OutputFolder', 'week2_1s');
%     S = process_LocomToSleep(timeInterval=1, folderPath='week2_1s', selectAll=true, startZT=...);
%
%   Golden test: test_RawToActivity1s (checks both modes against the Python reference, beam_hmm/labconv.py).
%   Written for MATLAB R2025a+; also runs in Octave (no arguments block).
%
%   TYPICAL SEQUENCE (functions in beam_hmm/matlab: cd there or addpath it; one weekly raw file at a time)
%     1  After loading the tubes, record ~2 h, then screen the tubes and re-seat the flagged ones:
%          R = LoadingCheck('raw_positions_all64_week2.txt', 'Channels', 17:32);
%        After re-seating (here at ~2.3 h), record >= 2 h more and confirm (RESOLVED / STILL FLAGGED / NEW FLAG):
%          R2 = LoadingCheck('raw_positions_all64_week2.txt', 'Channels', 17:32, 'StartHour', 2.5, 'Previous', R);
%     2  After the recording, look at the whole file (QC table, channel status, four figures):
%          F = BeamFingerprint('raw_positions_all64_week2.txt', 'Channels', 17:32, 'SaveFolder', 'fp_week2');
%     3  Convert to activity per second (repairs + QC), one .txt per channel:
%          A = RawToActivity1s('raw_positions_all64_week2.txt', 'Channels', 17:32, 'OutputFolder', 'week2_1s');
%     4  Score sleep with the lab pipeline, unchanged (one experiment per folder; set startZT -- nothing checks it):
%          S = process_LocomToSleep(timeInterval=1, folderPath='week2_1s', selectAll=true, startZT=0);
%     Optional: EstimatePhaseStamp to check the phase (startZT) from behaviour; chi2p for the period in DD.
%     Once per installation: test_BeamFingerprint and test_RawToActivity1s on the week-1 raw file (golden tests).
%     (Step 4's sleeptrace + ironout score this output identically to beam_hmm, verified second by second.)
%   NOTE      the startZT you give process_LocomToSleep is NOT checked by anything; in the 2026 pilot several
%             weekly files started 30-40 min away from lights-on (computer and incubator clocks drifted).
%   Definitions of every term and flag (hub, spike, FROZEN, P1-F death rule, S/Q/W ...): beam_hmm/GLOSSARY.md
%   See also LOADINGCHECK, BEAMFINGERPRINT, PROCESS_LOCOMTOSLEEP, SLEEPTRACE, IRONOUT, PREPROCESS,
%            MULTIBEAM_LOCOMOTION.

p = inputParser;
p.addRequired('filePath', @(s) ischar(s) || isstring(s));
p.addParameter('Mode', 'beam_hmm', @(s) any(strcmpi(char(s), {'beam_hmm', 'legacy'})));
p.addParameter('Repair', true);
p.addParameter('KeepBeyondBeam1', true);
p.addParameter('Channels', 1:64);
p.addParameter('OutputFolder', '');
p.addParameter('FilePrefix', '');
p.addParameter('Overwrite', false);
p.addParameter('AtZeroFolder', '');
p.parse(filePath, varargin{:});
o = p.Results; filePath = char(filePath); mode = lower(char(o.Mode));
ch = o.Channels(:)';

%% read every data line (the -1 lines included; legacy mode needs them)
fid = fopen(filePath, 'r');
if fid < 0, error('RawToActivity1s:NoFile', 'Cannot open %s', filePath); end
hdr = cell(3, 1);
for i = 1:3, hdr{i} = strtrim(fgetl(fid)); end
C = textscan(fid, [repmat('%d8', 1, 64) '%f'], 'Delimiter', '\t', 'EmptyValue', NaN, 'CollectOutput', true);
fclose(fid);
A = C{1}; tAll = C{2};
timed = ~isnan(tAll);
if any(any(A(~timed, :) ~= -1))
    error('RawToActivity1s:Format', 'Untimed lines contain values other than -1: the file format has changed.');
end
t = tAll(timed); stepS = median(diff(t));
info = struct('File', filePath, 'Header', {hdr}, 'Lines', size(A, 1), 'TimedLines', numel(t), 'StepSeconds', stepS, ...
              'IrregularSteps', sum(abs(diff(t) - stepS) > 1e-6));

%% per channel
nC = numel(ch); counts = cell(1, nC); zeros1s = cell(1, nC);
[hubV, hubShare, hubRet, hubDead, hubN, spikes, rangeB, ch24, at0, ent0] = deal(nan(nC, 1));
nSec = floor(t(end)) + 1; secAll = floor(t) + 1; nPerSec = accumarray(secAll, 1, [nSec 1]);
days = (t(end) - t(1)) / 86400;
flags = repmat({''}, nC, 1);
for k = 1:nC
    x = double(A(timed, ch(k)));
    % ---- QC on the timed readings (both modes)
    xs = x; nsp = 0;
    if o.Repair, [xs, nsp] = spikeFilter(x); end
    [xr, hi] = hubRepair(xs);
    hubV(k) = hi.hub; hubShare(k) = hi.share; hubRet(k) = hi.returns; hubDead(k) = hi.replaced > 0; hubN(k) = hi.replaced;
    spikes(k) = nsp; rangeB(k) = max(xr) - min(xr);
    ch24(k) = sum(diff(xr(t - t(1) < 86400)) ~= 0); at0(k) = mean(x == 0);
    f = {};
    if ch24(k) < 100, f{end + 1} = sprintf('EMPTY/DEAD: %d position changes in the first 24 h', ch24(k)); end
    if hubDead(k)
        if strcmp(mode, 'beam_hmm') && o.Repair
            f{end + 1} = sprintf('HUB %g declared dead (%.1f%% of readings replaced)', hubV(k), 100 * hubN(k) / numel(x));
        else
            f{end + 1} = sprintf('HUB %g DETECTED, NOT repaired in this mode (%.1f%% of readings): its flicker counts as movement', hubV(k), 100 * hubN(k) / numel(x));
        end
        if any(hubV(k) == [1 2 15 16])
            f{end + 1} = 'END-BEAM HUB: position at that end unknowable (food or debris blocking the end beam?)';
        end
    end
    if ch24(k) >= 100 && rangeB(k) < 8, f{end + 1} = sprintf('NARROW RANGE: %g beams spanned', rangeB(k)); end
    % ---- position 0 (feeding proxy in a correctly oriented tube): same repairs as the activity
    if o.Repair, xz = xr; else, xz = x; end
    z = xz == 0; ent0(k) = sum(z(2:end) & ~z(1:end - 1)) / days;
    zs = accumarray(secAll, double(z), [nSec 1]) ./ nPerSec; zs(nPerSec == 0) = NaN; zeros1s{k} = zs;
    flags{k} = strjoin(f, '; ');
    % ---- activity
    if strcmp(mode, 'legacy')
        a = legacyPreProcess(double(A(:, ch(k))));
        d = [abs(diff(a)); 0]; d(end) = d(end - 1);
        n5 = floor(numel(d) / 5) * 5;
        counts{k} = sum(reshape(d(1:n5), 5, []), 1)';
    else
        if o.Repair, xa = xr; else, xa = x; end
        if ~o.KeepBeyondBeam1
            for i = 2:numel(xa)
                if xa(i) == 0, xa(i) = xa(i - 1); end
            end
        end
        d = [0; abs(diff(xa))]; sec = floor(t);
        n = floor(t(end)) + 1;
        counts{k} = accumarray(sec + 1, d, [n 1]);
    end
end
nS = min(cellfun(@numel, counts));
Counts = zeros(nS, nC); AtZero = nan(nS, nC);
for k = 1:nC, Counts(:, k) = counts{k}(1:nS); m = min(nS, nSec); AtZero(1:m, k) = zeros1s{k}(1:m); end

%% QC printout and table
status = repmat({'INCLUDED'}, nC, 1); status(ch24 < 100) = {'EMPTY/DEAD'};
fprintf('RawToActivity1s: %s | mode %s | repair %d | keep beyond beam 1 %d | %d timed lines, step %.3f s\n', ...
        filePath, mode, o.Repair && strcmp(mode, 'beam_hmm'), o.KeepBeyondBeam1, numel(t), stepS);
for k = 1:nC
    if ~isempty(flags{k}), fprintf('   ch%-3d %s\n', ch(k), flags{k}); end
end
names = {'Channel', 'Status', 'Flags', 'HubValue', 'HubShare', 'HubReturns', 'HubDeclaredDead', 'HubReadingsReplaced', ...
         'SpikesFixed', 'BeamsSpanned', 'Changes24h', 'ShareAtZero', 'EntriesToZeroPerDay'};
cols = {ch', status, flags, hubV, hubShare, hubRet, logical(hubDead), hubN, spikes, rangeB, ch24, at0, ent0};
if exist('table', 'file') || exist('table', 'builtin')
    QC = table(cols{:}, 'VariableNames', names);
else
    QC = cell2struct(cols, names, 2);
end

%% optional output files
if ~isempty(o.OutputFolder)
    if ~exist(o.OutputFolder, 'dir'), mkdir(o.OutputFolder); end
    pre = o.FilePrefix;
    if isempty(pre), [~, nm] = fileparts(filePath); pre = [nm '_ch']; end
    for k = 1:nC
        fn = fullfile(o.OutputFolder, sprintf('%s%d.txt', pre, ch(k)));
        if exist(fn, 'file') && ~o.Overwrite
            error('RawToActivity1s:Exists', '%s exists; refusing to append or overwrite (set ''Overwrite'', true).', fn);
        end
        fid = fopen(fn, 'w'); fprintf(fid, '%d \r\n', Counts(:, k)); fclose(fid);
    end
end
if ~isempty(o.AtZeroFolder)
    if ~isempty(o.OutputFolder) && strcmp(fullfile(o.AtZeroFolder, '.'), fullfile(o.OutputFolder, '.'))
        error('RawToActivity1s:AtZeroFolder', 'AtZeroFolder must differ from OutputFolder (process_LocomToSleep reads every file there).');
    end
    if ~exist(o.AtZeroFolder, 'dir'), mkdir(o.AtZeroFolder); end
    pre = o.FilePrefix;
    if isempty(pre), [~, nm] = fileparts(filePath); pre = [nm '_ch']; end
    for k = 1:nC
        fn = fullfile(o.AtZeroFolder, sprintf('%s%d_at0.txt', pre, ch(k)));
        if exist(fn, 'file') && ~o.Overwrite
            error('RawToActivity1s:Exists', '%s exists; refusing to overwrite (set ''Overwrite'', true).', fn);
        end
        fid = fopen(fn, 'w'); fprintf(fid, '%.4f\r\n', AtZero(:, k)); fclose(fid);
    end
end
out = struct('Counts', Counts, 'AtZero', AtZero, 'Channels', ch, 'Mode', mode, 'Info', info, 'QC', QC, ...
             'Settings', struct('Repair', o.Repair, 'KeepBeyondBeam1', o.KeepBeyondBeam1));
end

%% ======================================================================== local functions
function a = legacyPreProcess(a)
% PreProcess.m on one column (all lines). The row bound is checked first (PreProcess checked it after
% indexing, which errors on a column of only 0s and 1s); otherwise identical.
n = numel(a); i = 1;
while i <= n && (a(i) == 0 || a(i) == 1), i = i + 1; end
if i <= n, a(1:i) = a(i); end
for k = 2:n
    if a(k) == 0 || a(k) == -1 || (a(k) == 1 && a(k - 1) ~= 2)
        a(k) = a(k - 1);
    end
end
end

function [x, n] = spikeFilter(x)
a = x(1:end - 2); b = x(2:end - 1); c = x(3:end);
idx = find(a == c & abs(b - a) >= 4) + 1;
x(idx) = x(idx - 1); n = numel(idx);
end

function [x, info] = hubRepair(x)
info = struct('hub', NaN, 'share', 0, 'returns', 0, 'replaced', 0);
big = find(abs(diff(x)) >= 4) + 1;
if isempty(big), return; end
[v, ~, j] = unique(x(big)); cnt = accumarray(j, 1);
[nm, im] = max(cnt); hub = v(im);
info.hub = hub; info.share = nm / numel(big);
starts = big(x(big) == hub); ret = false(numel(starts), 1);
for s = 1:numel(starts)
    i = starts(s); k = i;
    while k + 1 <= numel(x) && x(k + 1) == hub, k = k + 1; end
    ret(s) = (k + 1 <= numel(x)) && x(k + 1) == x(i - 1);
end
info.returns = mean(ret);
if info.share >= 0.30 && info.returns >= 0.20
    at = x == hub; info.replaced = sum(at);
    if at(1)
        nz = find(~at, 1);
        if ~isempty(nz), x(1) = x(nz); end
        at(1) = false;
    end
    idx = (1:numel(x))'; idx(at) = 0; idx = cummax(idx); x = x(idx);
end
end
