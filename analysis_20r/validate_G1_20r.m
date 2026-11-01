function V = validate_G1_20r(goldenFile, codeDir, outFile)
% VALIDATE_G1_20R  Run the hmm_audit golden bridging cases (G1_sleep_bouts.json,
% made by sqw_ref.sleep_bouts in Python) through the MATLAB pipeline's own
% scoring chain -- sleeptrace -> ironout -> segment_sleep_phases, as
% process_LocomToSleep calls them -- and compare sleep bouts case by case.
%
% Golden starts are 0-based; they are shifted by +1 here. Each case is
% reported twice for MATLAB: the bouts segment_sleep_phases keeps (it drops a
% bout touching the first or last sample as censored), and every run of
% sleep in the scored wave before that censoring, so a disagreement can be
% placed on bridging or on edge handling. Nothing is changed on either side.

addpath(codeDir, '-begin');
cleanupPath = onCleanup(@() rmpath(codeDir));
for f = {'sleeptrace','ironout','segment_sleep_phases'}
    assert(isequal(fileparts(which(f{1})), codeDir), '%s resolves to %s', f{1}, which(f{1}));
end

G = jsondecode(fileread(goldenFile));
names = fieldnames(G);
L = {sprintf('G1 golden validation, %s; code %s', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm')), codeDir)};
V = struct('case', {}, 'agreeBridged', {}, 'agreeUnbridged', {});
for i = 1:numel(names)
    c = G.(names{i});
    x = double(c.x(:));
    aw = sleeptrace(x, 1);
    for bridged = [true false]
        if bridged
            w = ironout(aw, 1, 1); ref = c.bouts_bridged; lab = 'bridged';
        else
            w = aw; ref = c.bouts_unbridged; lab = 'unbridged';
        end
        ref = reshape(ref, [], 2); ref(:, 1) = ref(:, 1) + 1;   % 0-based -> 1-based
        seg = segment_sleep_phases(w, timeInterval=1);
        sb  = seg.stats.sleep_bouts;
        kept = [sb.starts(:), sb.durations(:)];
        all0 = runsOf(w == 0);
        agree = isequal(sortrows(kept), sortrows(ref));
        L{end+1} = sprintf('%-26s %-9s python %-22s matlab kept %-22s matlab before censoring %-22s %s', ...
            names{i}, lab, fmt(ref), fmt(kept), fmt(all0), ternary(agree, 'AGREE', 'DIFFER')); %#ok<AGROW>
        if bridged, V(i).agreeBridged = agree; else, V(i).agreeUnbridged = agree; end
    end
    V(i).case = names{i};
end
txt = strjoin(L, newline);
fprintf('%s\n', txt);
fid = fopen(outFile, 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
end

function r = runsOf(m)
d = diff([false; m(:); false]);
s = find(d == 1); e = find(d == -1) - 1;
r = [s, e - s + 1];
r = r(r(:, 2) >= 300, :);                   % sleep: >= 300 s, as sleeptrace scores it
end

function s = fmt(b)
if isempty(b), s = '[]'; return, end
s = strjoin(arrayfun(@(k) sprintf('[%d %d]', b(k,1), b(k,2)), 1:size(b,1), 'UniformOutput', false), ' ');
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
