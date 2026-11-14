function R = LoadingCheck(filePath, varargin)
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
%     LOW ACTIVITY      < 20 position changes: a dead or stuck fly, or a blocked tube.
%   A fly resting on one beam is NOT flagged by "reads one beam a lot"; only the jump signature counts.
%   These are the beam_hmm hub criteria (PRESPEC P1-C) with a minimum jump count for short recordings.
%
%   NAME-VALUE: 'Hours' (default 2), 'Channels' (default 1:64).
%   Python twin: beam_hmm/labconv.py loading_check (validated against the full-week QC of a 9-week run).

p = inputParser;
p.addRequired('filePath', @(s) ischar(s) || isstring(s));
p.addParameter('Hours', 2); p.addParameter('Channels', 1:64);
p.parse(filePath, varargin{:}); o = p.Results; filePath = char(filePath);
fid = fopen(filePath, 'r');
if fid < 0, error('LoadingCheck:NoFile', 'Cannot open %s', filePath); end
for i = 1:3, fgetl(fid); end
nLines = ceil(o.Hours * 3600 / 0.2) + 1000;            % enough lines for the window at <= 5 Hz, -1 lines included
C = textscan(fid, [repmat('%d8', 1, 64) '%f'], nLines, 'Delimiter', '\t', 'EmptyValue', NaN, 'CollectOutput', true);
fclose(fid);
A = C{1}; t = C{2}; keep = ~isnan(t); A = A(keep, :); t = t(keep);
m = (t - t(1)) < o.Hours * 3600; A = double(A(m, :));
if (t(find(m, 1, 'last')) - t(1)) < 0.95 * o.Hours * 3600
    warning('LoadingCheck:Short', 'File holds only %.2f h; checked what there is.', (t(find(m, 1, 'last')) - t(1)) / 3600);
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
    if nchg(k) < 20, f{end + 1} = 'LOW ACTIVITY'; end
    flags{k} = strjoin(f, '; ');
end
bad = find(~cellfun(@isempty, flags));
fprintf('LoadingCheck: %s, first %.1f h, %d channels\n', filePath, o.Hours, nC);
if isempty(bad)
    fprintf('  no tube flagged\n');
else
    fprintf('  RE-SEAT / CHECK %d tube(s):\n', numel(bad));
    for k = bad(:)', fprintf('    ch%-3d %s\n', ch(k), flags{k}); end
end
names = {'Channel', 'Flags', 'Changes', 'BigJumps', 'HubValue', 'HubShare', 'HubReturns', 'BeamsSpanned'};
cols = {ch', flags, nchg, nbig, hubV, hubSh, hubRet, span};
if exist('table', 'file') || exist('table', 'builtin'), R = table(cols{:}, 'VariableNames', names);
else, R = cell2struct(cols, names, 2); end
end
