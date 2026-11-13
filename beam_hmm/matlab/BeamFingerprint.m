function out = BeamFingerprint(filePath, varargin)
%BEAMFINGERPRINT  QC table and four "fingerprint" figures for one raw multibeam file.
%
%   OUT = BEAMFINGERPRINT(FILEPATH) reads one raw multibeam position file (one column per tube,
%   one row per sample, a time column last), applies the beam_hmm preprocessing, and draws:
%     1. activity heatmap: channels x hour of the file (share of samples with a position change)
%     2. raw position traces of selected channels over one day, dark phase shaded
%     3. 10-minute zoom of the same channels
%     4. occupancy by ZT hour: share of time at each beam, per channel (small multiples)
%   It returns a per-channel QC table (out.QC) and the figure handles (out.figs).
%
%   PREPROCESSING (identical to beam_hmm/pipeline.py, PRESPEC P1 and P1-C; see beam_hmm/PRESPEC.md)
%     - lines whose 64 position fields are all -1 carry no time stamp: they are DROPPED, not filled
%       (the true sampling step is then 0.4 s although the header says 5 Hz; the step is measured,
%       never assumed: out.Info.StepSeconds)
%     - value 0 is kept as its own position, "beyond beam 1" (it is always entered from beam 1)
%     - spike filter, one pass: a sample B with equal neighbours A (A-B-A) and |B-A| >= 4 is set to A
%     - dead-hub repair: the value receiving the most >= 4-beam jumps is a HUB if it receives >= 30% of
%       them AND >= 20% of the runs entered that way return to the pre-jump value; then EVERY sample
%       reading the hub is replaced by the last non-hub value (the hub beam is treated as dead, as the
%       hardware already treats beam 15)
%     - a channel is LIVE if it has >= 100 position changes in its first 24 h, and ARTIFACT-HEAVY if
%       > 5% of its changes are >= 4-beam jumps after both repairs (P1-D)
%   Thresholds are those fixed in beam_hmm/PRESPEC.md; change them only with a new entry there.
%
%   THIS IS A DESCRIPTIVE TOOL. It computes no test. Its figures are what a session looks at FIRST, so
%   anything seen in them is post hoc for any later analysis of the same file.
%
%   NAME-VALUE OPTIONS
%     'Channels'        channels to analyse (default 1:64; dead/empty ones are reported, not plotted)
%     'LightsOnHour'    hours into the file at which lights come on (default 0: file starts at lights-on)
%     'Photoperiod'     [light dark] hours (default [12 12])
%     'TraceDay'        day of the file for figure 2 (1-based; default 2)
%     'TraceChannels'   channels for figures 2-3 (default: up to 8 live, non-artifact channels)
%     'ZoomStartHours'  start of the 10-min zoom, hours into TraceDay (default 3)
%     'ZoomMinutes'     zoom length (default 10)
%     'OccupancyHours'  [start end] hours into the file used for figure 4 (default: whole file)
%     'OccupancyChannels' channels for figure 4 (default: live, non-artifact channels)
%     'Repair'          apply spike filter and hub repair (default true)
%     'SaveFolder'      if non-empty, each figure is also saved there as PNG (default '')
%     'Visible'         figure visibility 'on'/'off' (default 'on')
%     'Plot'            false = QC table and status printout only, no figures (default true)
%
%   EXAMPLE
%     out = BeamFingerprint('raw_positions_all64_week2.txt', 'Channels', 17:32, 'SaveFolder', 'fp_week2');
%     disp(out.QC)
%
%   Written for MATLAB R2025a+; also runs in Octave (no arguments block, no tiledlayout). The figure
%   theme is pinned to light (CLAUDE.md: an unpinned R2026a figure can come out dark).

p = inputParser;
p.addRequired('filePath', @(s) ischar(s) || isstring(s));
p.addParameter('Channels', 1:64);
p.addParameter('LightsOnHour', 0);
p.addParameter('Photoperiod', [12 12]);
p.addParameter('TraceDay', 2);
p.addParameter('TraceChannels', []);
p.addParameter('ZoomStartHours', 3);
p.addParameter('ZoomMinutes', 10);
p.addParameter('OccupancyHours', []);
p.addParameter('OccupancyChannels', []);
p.addParameter('Repair', true);
p.addParameter('SaveFolder', '');
p.addParameter('Visible', 'on');
p.addParameter('Plot', true);
p.parse(filePath, varargin{:});
o = p.Results;
filePath = char(filePath);

%% 1. read ---------------------------------------------------------------
[X, t, hdr] = readBeamFile(filePath);
ch = o.Channels(:)';
X = X(:, ch);
dt = diff(t);
stepS = median(dt);
info = struct('File', filePath, 'Header', {hdr}, 'NSamples', numel(t), 'StepSeconds', stepS, ...
              'IrregularSteps', sum(abs(dt - stepS) > 1e-6), 'DurationHours', (t(end) - t(1)) / 3600);
fprintf('BeamFingerprint: %s\n  %d timed samples, step %.3f s (%d irregular), %.1f h, %d channels\n', ...
        filePath, info.NSamples, stepS, info.IrregularSteps, info.DurationHours, numel(ch));

%% 2. per-channel preprocessing and QC --------------------------------------
nC = numel(ch);
P = X;
[nVals, vMin, vMax, nChg, nChg24, fBigRaw, fBigFinal, spikes, hub, hubShare, hubRet, hubN, ...
 silH, silStartH] = deal(nan(nC, 1));
hubOn = false(nC, 1);
first24 = t - t(1) < 86400;
for k = 1:nC
    x = double(X(:, k));
    nVals(k) = numel(unique(x)); vMin(k) = min(x); vMax(k) = max(x);
    fBigRaw(k) = bigFrac(x);
    if o.Repair
        [x, spikes(k)] = spikeFilter(x);
        [x, hi] = hubRepair(x);
        hub(k) = hi.hub; hubShare(k) = hi.share; hubRet(k) = hi.returns; hubN(k) = hi.replaced;
        hubOn(k) = hi.replaced > 0;
    end
    P(:, k) = x;
    c = [false; diff(x) ~= 0];
    nChg(k) = sum(c); nChg24(k) = sum(c & first24);
    fBigFinal(k) = bigFrac(x);
    tc = t(c);
    edges = [t(1); tc; t(end)];
    [g, gi] = max(diff(edges));
    silH(k) = g / 3600; silStartH(k) = (edges(gi) - t(1)) / 3600;
end
live = nChg24 >= 100;
artifact = live & fBigFinal > 0.05;
status = repmat({'included'}, nC, 1);
status(~live) = {'empty/dead at load'};
status(artifact) = {'artifact-heavy'};
good = live & ~artifact;
% per-channel flags, written out so that nothing about a channel has to be inferred from the figures
flags = repmat({''}, nC, 1);
for k = 1:nC
    f = {};
    if ~live(k)
        if nChg(k) == 0
            f{end + 1} = sprintf('EMPTY: no position change in the whole file (constant %g)', vMin(k));
        else
            f{end + 1} = sprintf('EMPTY/DEAD: %d changes in the first 24 h (< 100)', nChg24(k));
        end
    end
    if live(k) && (vMax(k) - vMin(k) < 8)
        f{end + 1} = sprintf('NARROW RANGE: only values %g-%g (%d distinct): check the hardware before use', vMin(k), vMax(k), nVals(k));
    end
    if hubOn(k)
        f{end + 1} = sprintf('HUB %g declared dead (%.1f%% of samples replaced)', hub(k), 100 * hubN(k) / numel(t));
    end
    if artifact(k)
        f{end + 1} = sprintf('ARTIFACT-HEAVY: %.1f%% of changes are >= 4-beam jumps after repair', 100 * fBigFinal(k));
    end
    if live(k) && silH(k) >= 12
        f{end + 1} = sprintf('SILENT %.1f h from hour %.1f (death or handling? check before cutting)', silH(k), silStartH(k));
    end
    flags{k} = strjoin(f, '; ');
end
fprintf('\n  CHANNEL STATUS (%d channels)\n', nC);
groups = {'EMPTY / DEAD AT LOAD', ~live; 'ARTIFACT-HEAVY (reported, not pooled)', artifact; 'INCLUDED', good};
for gI = 1:size(groups, 1)
    idx = find(groups{gI, 2});
    fprintf('  %s: %d', groups{gI, 1}, numel(idx));
    if isempty(idx), fprintf('\n'); continue; end
    fprintf('  -> channels %s\n', strjoin(arrayfun(@num2str, ch(idx), 'UniformOutput', false), ' '));
    for k = idx(:)'
        if ~isempty(flags{k}), fprintf('     ch%-3d %s\n', ch(k), flags{k}); end
    end
end
fprintf('\n');
qcNames = {'Channel', 'Status', 'Flags', 'NValues', 'Min', 'Max', 'Changes24h', 'Changes', ...
           'SpikesFixed', 'HubValue', 'HubShare', 'HubReturns', 'HubDeclaredDead', ...
           'HubSamplesReplaced', 'FracJumpsGE4_raw', 'FracJumpsGE4_final', ...
           'LongestSilence_h', 'LongestSilenceStart_h'};
qcCols = {ch', status, flags, nVals, vMin, vMax, nChg24, nChg, spikes, hub, hubShare, hubRet, hubOn, hubN, ...
          fBigRaw, fBigFinal, silH, silStartH};
if exist('table', 'file') || exist('table', 'builtin')
    QC = table(qcCols{:}, 'VariableNames', qcNames);
else                                   % Octave: struct of columns
    QC = cell2struct(qcCols, qcNames, 2);
end

%% 3. figures -----------------------------------------------------------------
zt = mod((t - t(1)) / 3600 - o.LightsOnHour, 24);
hr = floor((t - t(1)) / 3600);
if exist('gobjects', 'builtin') || exist('gobjects', 'file'), figs = gobjects(0); else, figs = []; end
liveIdx = find(live);

% 3.1 activity heatmap
if o.Plot && ~isempty(liveIdx)
    nH = max(hr) + 1;
    A = nan(numel(liveIdx), nH);
    n = accumarray(hr + 1, 1, [nH 1]);
    for j = 1:numel(liveIdx)
        c = [false; diff(P(:, liveIdx(j))) ~= 0];
        A(j, :) = accumarray(hr + 1, double(c), [nH 1]) ./ max(n, 1);
    end
    fh = newFigure(o.Visible, [100 100 1100 520], 'Activity heatmap');
    ax = axes(fh);
    imagesc(ax, 0:nH - 1, 1:numel(liveIdx), A);
    caxis(ax, [0 max(quantile99(A(:)), eps)]);
    colormap(ax, parula); cb = colorbar(ax); ylabel(cb, 'share of samples with a position change');
    hold(ax, 'on');
    for h = o.LightsOnHour:12:nH
        plot(ax, [h h] - 0.5, [0.5 numel(liveIdx) + 0.5], '-', 'Color', [1 1 1] * (mod((h - o.LightsOnHour) / 12, 2) == 1), 'LineWidth', 0.5);
    end
    set(ax, 'YTick', 1:numel(liveIdx), 'YTickLabel', arrayfun(@num2str, ch(liveIdx), 'UniformOutput', false), 'FontSize', 7);
    xlabel(ax, 'hour from file start (lines every 12 h at light transitions)'); ylabel(ax, 'channel');
    excl = ch(~live);
    ttl = sprintf('%s: activity by channel and hour (live channels only)', shortName(filePath));
    if ~isempty(excl)
        ttl = sprintf('%s\nnot shown, empty/dead: %s', ttl, strjoin(arrayfun(@num2str, excl, 'UniformOutput', false), ' '));
    end
    title(ax, ttl, 'Interpreter', 'none');
    figs(end + 1) = fh;
end

% channels for traces
tr = o.TraceChannels;
if isempty(tr)
    gi = find(good); tr = ch(gi(1:min(8, numel(gi))));
end
[~, trIdx] = ismember(tr, ch); trIdx = trIdx(trIdx > 0);

% 3.2 traces over one day
dayStart = (o.TraceDay - 1) * 86400 + o.LightsOnHour * 3600;
m = (t - t(1)) >= dayStart & (t - t(1)) < dayStart + 86400;
if o.Plot && any(m) && ~isempty(trIdx)
    fh = newFigure(o.Visible, [100 100 1100 140 * numel(trIdx) + 80], 'Raw traces, one day');
    for j = 1:numel(trIdx)
        ax = subplot(numel(trIdx), 1, j, 'Parent', fh); hold(ax, 'on');
        x0 = ((t(m) - t(1)) - dayStart) / 3600;
        patch(ax, [o.Photoperiod(1) 24 24 o.Photoperiod(1)], [-0.5 -0.5 16.5 16.5], [0.85 0.85 0.85], 'EdgeColor', 'none');
        stairs(ax, x0, P(m, trIdx(j)), 'k', 'LineWidth', 0.3);
        ylim(ax, [-0.5 16.5]); xlim(ax, [0 24]); ylabel(ax, sprintf('ch%d', ch(trIdx(j))));
        if j < numel(trIdx), set(ax, 'XTickLabel', []); end
    end
    xlabel(ax, sprintf('ZT (h), day %d of the file; grey = dark; position 0 = beyond beam 1', o.TraceDay));
    figs(end + 1) = fh;
end

% 3.3 zoom
z0 = dayStart + o.ZoomStartHours * 3600;
m = (t - t(1)) >= z0 & (t - t(1)) < z0 + o.ZoomMinutes * 60;
if o.Plot && any(m) && ~isempty(trIdx)
    nz = min(4, numel(trIdx));
    fh = newFigure(o.Visible, [100 100 1100 160 * nz + 80], 'Zoom');
    for j = 1:nz
        ax = subplot(nz, 1, j, 'Parent', fh);
        stairs(ax, (t(m) - t(1)) - z0, P(m, trIdx(j)), 'k', 'LineWidth', 0.6);
        ylim(ax, [-0.5 16.5]); ylabel(ax, sprintf('ch%d', ch(trIdx(j))));
        if j < nz, set(ax, 'XTickLabel', []); end
    end
    xlabel(ax, sprintf('s (day %d, ZT %g + %g min)', o.TraceDay, o.ZoomStartHours, o.ZoomMinutes));
    figs(end + 1) = fh;
end

% 3.4 occupancy by ZT hour
oc = o.OccupancyChannels;
if isempty(oc), oc = ch(good); end
[~, ocIdx] = ismember(oc, ch); ocIdx = ocIdx(ocIdx > 0);
hrsSel = true(size(t));
if ~isempty(o.OccupancyHours)
    hh = (t - t(1)) / 3600; hrsSel = hh >= o.OccupancyHours(1) & hh < o.OccupancyHours(2);
end
if o.Plot && ~isempty(ocIdx)
    zh = floor(zt(hrsSel));
    nCol = min(8, numel(ocIdx)); nRow = ceil(numel(ocIdx) / nCol);
    fh = newFigure(o.Visible, [100 100 min(1700, 210 * nCol + 120) 300 * nRow + 60], 'Occupancy by ZT');
    for j = 1:numel(ocIdx)
        x = P(hrsSel, ocIdx(j));
        M = accumarray([x + 1, zh + 1], 1, [17 24]);
        M = M ./ max(sum(M, 1), 1);
        ax = subplot(nRow, nCol, j, 'Parent', fh);
        imagesc(ax, 0:23, 0:16, M); set(ax, 'YDir', 'normal'); caxis(ax, [0 0.5]); colormap(ax, hot);
        hold(ax, 'on'); plot(ax, [11.5 11.5], [-0.5 16.5], 'c-', 'LineWidth', 0.8);
        title(ax, sprintf('ch%d', ch(ocIdx(j)))); xlabel(ax, 'ZT h');
        if mod(j - 1, nCol) == 0, ylabel(ax, 'beam (0 = beyond beam 1)'); end
    end
    if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
        sgtitle(fh, 'Occupancy by ZT hour (descriptive)');
    end
    figs(end + 1) = fh;
end

%% 4. save -----------------------------------------------------------------
if ~isempty(o.SaveFolder)
    if ~exist(o.SaveFolder, 'dir'), mkdir(o.SaveFolder); end
    names = {'activity_heatmap', 'traces_day', 'traces_zoom', 'occupancy_zt'};
    for j = 1:numel(figs)
        print(figs(j), fullfile(o.SaveFolder, sprintf('%s_%s.png', shortName(filePath), names{j})), '-dpng', '-r150');
    end
    if isTableObj(QC)
        writetable(QC, fullfile(o.SaveFolder, sprintf('%s_QC.csv', shortName(filePath))));
    end
end
out = struct('QC', QC, 'Info', info, 'figs', figs, 'Positions', P, 'Time', t, 'Channels', ch);
end

%% ===================================================================== local functions

function [X, t, hdr] = readBeamFile(f)
% Header: 3 lines. Data: 64 tab-separated integers, then a time field (empty on the all -1 lines).
fid = fopen(f, 'r');
if fid < 0, error('BeamFingerprint:NoFile', 'Cannot open %s', f); end
hdr = cell(3, 1);
for i = 1:3, hdr{i} = strtrim(fgetl(fid)); end
fmt = [repmat('%d8', 1, 64) '%f'];          % positions fit int8 (-1..16); keeps memory ~200 MB/week
C = textscan(fid, fmt, 'Delimiter', '\t', 'EmptyValue', NaN, 'CollectOutput', true, 'ReturnOnError', false);
fclose(fid);
A = C{1}; t = C{2};
keep = ~isnan(t);
if any(any(A(~keep, :) ~= -1))
    error('BeamFingerprint:Format', 'Untimed lines contain values other than -1: the format has changed.');
end
X = int16(A(keep, :));
t = t(keep);
end

function f = bigFrac(x)
d = abs(diff(x)); c = d ~= 0;
if any(c), f = mean(d(c) >= 4); else, f = NaN; end
end

function [x, n] = spikeFilter(x)
a = x(1:end - 2); b = x(2:end - 1); c = x(3:end);
idx = find(a == c & abs(b - a) >= 4) + 1;
x(idx) = x(idx - 1);
n = numel(idx);
end

function [x, info] = hubRepair(x)
info = struct('hub', NaN, 'share', 0, 'returns', 0, 'replaced', 0);
d = diff(x); big = find(abs(d) >= 4) + 1;
if isempty(big), return; end
[v, ~, j] = unique(x(big)); n = accumarray(j, 1);
[nm, im] = max(n); hub = v(im);
info.hub = hub; info.share = nm / numel(big);
starts = big(x(big) == hub);
ret = false(numel(starts), 1);
for s = 1:numel(starts)
    i = starts(s); k = i;
    while k + 1 <= numel(x) && x(k + 1) == hub, k = k + 1; end
    ret(s) = (k + 1 <= numel(x)) && x(k + 1) == x(i - 1);
end
info.returns = mean(ret);
if info.share >= 0.30 && info.returns >= 0.20
    at = x == hub;
    info.replaced = sum(at);
    if at(1)
        nz = find(~at, 1);
        if ~isempty(nz), x(1) = x(nz); end
        at(1) = false;
    end
    idx = (1:numel(x))'; idx(at) = 0;
    idx = cummax(idx);
    x = x(idx);
end
end

function q = quantile99(v)
v = sort(v(~isnan(v)));
if isempty(v), q = 0; return; end
q = v(max(1, ceil(0.99 * numel(v))));
end

function tf = isTableObj(x)
tf = isa(x, 'table');
end

function fh = newFigure(vis, pos, name)
fh = figure('Color', 'w', 'Position', pos, 'Name', name, 'Visible', vis);
try
    theme(fh, 'light');      % before any axes, so children inherit (R2025a+; absent in Octave)
catch
end
end

function s = shortName(f)
[~, s] = fileparts(char(f));
end
