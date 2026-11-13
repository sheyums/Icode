D = 'C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data\per0 dd';
T = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\per0dd_stage';
if isfolder(T), error('stage folder exists'); end
mkdir(T);
all = dir(fullfile(D, '*.txt'));
names = string({all.name});
cut = names(endsWith(names, '_from1stLD.txt'));
twins = erase(cut, '_from1stLD.txt');
plain = names(~endsWith(names, '_from1stLD.txt') & ~ismember(names, twins));
use = [cut plain];
fprintf('cut %d, plain kept %d, plain dropped (twins) %d, total %d\n', numel(cut), numel(plain), numel(twins(ismember(twins, names))), numel(use));
md = java.security.MessageDigest.getInstance('SHA-256');
fp = strings(numel(use), 1); n = zeros(numel(use), 1); seqHash = strings(numel(use), 1);
for i = 1:numel(use)
    src = fullfile(D, use(i));
    copyfile(src, fullfile(T, use(i)));
    fid = fopen(src); b = fread(fid, inf, '*uint8'); fclose(fid);
    md.reset(); h = typecast(md.digest(b), 'uint8'); fp(i) = string(upper(reshape(dec2hex(h(1:8))', 1, [])));
    x = readmatrix(src, 'FileType', 'text'); x = x(:, 1); x = x(~isnan(x));
    n(i) = numel(x);
    md.reset(); h = typecast(md.digest(typecast(int32(x), 'uint8')), 'uint8');
    seqHash(i) = string(reshape(dec2hex(h(1:8))', 1, []));
end
[u, ~, j] = unique(seqHash); c = accumarray(j, 1);
fprintf('distinct recordings by content: %d of %d\n', numel(u), numel(use));
% also check every recording's first day against every other (catches overlapping cuts)
firstDay = cell(numel(use), 1);
for i = 1:numel(use)
    x = readmatrix(fullfile(D, use(i)), 'FileType', 'text'); x = x(:, 1); x = x(~isnan(x));
    firstDay{i} = char(x(1:86400).' + 48); full{i} = char(x.' + 48); %#ok<SAGROW>
end
hits = 0;
for i = 1:numel(use)
    for k = 1:numel(use)
        if i ~= k && ~isempty(strfind(full{k}, firstDay{i}))
            fprintf('  first day of %s found inside %s\n', use(i), use(k)); hits = hits + 1;
        end
    end
end
fprintf('cross-file first-day matches: %d\n', hits);
tbl = table(use(:), n, n/86400, fp, 'VariableNames', {'File','Samples','Days','SHA256prefix'});
writetable(tbl, fullfile(fileparts(T), 'per0dd_manifest.csv'));
disp(tbl)
