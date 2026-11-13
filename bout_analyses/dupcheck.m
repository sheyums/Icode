D = 'C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data\per0 dd';
f = dir(fullfile(D, '*_from1stLD.txt'));
for i = 1:numel(f)
    cutName = f(i).name;
    plainName = erase(cutName, '_from1stLD.txt');
    c = readmatrix(fullfile(D, cutName), 'FileType', 'text');
    c = c(:, 1); c = c(~isnan(c));
    if ~isfile(fullfile(D, plainName))
        fprintf('%-45s  n=%7d (%.2f d)  no plain twin\n', cutName, numel(c), numel(c)/86400);
        continue
    end
    p = readmatrix(fullfile(D, plainName), 'FileType', 'text');
    p = p(:, 1); p = p(~isnan(p));
    k = strfind(char(p.' + 48), char(c.' + 48));
    fprintf('%-45s  cut n=%7d (%.2f d), plain n=%7d (%.2f d), cut found in plain at sample %s\n', ...
        cutName, numel(c), numel(c)/86400, numel(p), numel(p)/86400, mat2str(k));
end
