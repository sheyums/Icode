function T = align_cut_20r(cutDir, rawDir, outFile)
% ALIGN_CUT_20R  Locate each cut 1-s count file inside its raw *_from1stLD.txt
% source, to document the cut (offset from the start, length kept, samples left
% after it) so the same cut can be applied to another genotype.
%
%   cutDir  : folder of *_from1stLD.txtcutdata_*.txt files
%   rawDir  : folder of the matching raw *_from1stLD.txt files
%   outFile : text report to write
%
% A cut file matches at offset o when raw(o+1 : o+numel(cut)) equals cut
% exactly. Every exact match is listed, so an ambiguous placement (a long
% all-zero stretch, say) shows up rather than being resolved silently.

d = dir(fullfile(cutDir, '*_from1stLD.txtcutdata_*.txt'));
n = numel(d);
file = strings(n, 1); nRaw = nan(n, 1); nCut = nan(n, 1);
offset = nan(n, 1); after = nan(n, 1); nMatch = zeros(n, 1);
for i = 1:n
    cutName = d(i).name;
    rawName = regexprep(cutName, 'cutdata_.*\.txt$', '');
    file(i) = string(rawName);
    c = readmatrix(fullfile(cutDir, cutName));  c = c(:);
    r = readmatrix(fullfile(rawDir, rawName));  r = r(:);
    nCut(i) = numel(c); nRaw(i) = numel(r);
    m = numel(c);
    hits = [];
    if m <= numel(r)
        % candidate offsets from the first 600 samples, then confirm in full
        k = min(600, m);
        cand = strfind(r(1:end-m+k)', c(1:k)') - 1;
        for o = cand
            if isequal(r(o+1:o+m), c), hits(end+1) = o; end %#ok<AGROW>
        end
    end
    nMatch(i) = numel(hits);
    if ~isempty(hits)
        offset(i) = hits(1);
        after(i) = numel(r) - hits(1) - m;
    end
end

T = table(file, nRaw, nCut, offset, after, nMatch, ...
    'VariableNames', {'file','rawSamples','cutSamples','offsetSamples','samplesAfter','exactMatches'});
L = {sprintf('per0 DD cut alignment, %s', char(datetime('now','Format','yyyy-MM-dd HH:mm'))), ...
     sprintf('cut: %s', cutDir), sprintf('raw: %s', rawDir), ...
     sprintf('%-50s %9s %9s %9s %9s %7s  %s', 'raw file', 'raw', 'cut', 'offset', 'after', 'matches', 'offset/after (h)')};
for i = 1:n
    L{end+1} = sprintf('%-50s %9d %9d %9d %9d %7d  %.4f / %.4f', file(i), nRaw(i), nCut(i), ...
        offset(i), after(i), nMatch(i), offset(i)/3600, after(i)/3600); %#ok<AGROW>
end
txt = strjoin(L, newline);
fprintf('%s\n', txt);
fid = fopen(outFile, 'w'); fprintf(fid, '%s\n', txt); fclose(fid);
end
