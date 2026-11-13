D = 'C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data\CS ld';
T = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\csld_stage';
if isfolder(T), error('stage folder exists'); end
mkdir(T);
f = dir(fullfile(D, '*.txt')); use = string({f.name});
md = java.security.MessageDigest.getInstance('SHA-256');
n = zeros(numel(use),1); fp = strings(numel(use),1); sh = fp; X = cell(numel(use),1);
for i = 1:numel(use)
    src = fullfile(D, use(i)); copyfile(src, fullfile(T, use(i)));
    fid = fopen(src); b = fread(fid, inf, '*uint8'); fclose(fid);
    md.reset(); h = typecast(md.digest(b), 'uint8'); fp(i) = string(reshape(dec2hex(h(1:8))', 1, []));
    x = readmatrix(src, 'FileType', 'text'); x = x(:,1); x = x(~isnan(x)); X{i} = x; n(i) = numel(x);
    md.reset(); h = typecast(md.digest(typecast(int32(x), 'uint8')), 'uint8'); sh(i) = string(reshape(dec2hex(h(1:8))', 1, []));
end
fprintf('files %d, distinct contents %d\n', numel(use), numel(unique(sh)));
hits = 0;
for i = 1:numel(use), a = char(X{i}(1:86400).' + 48);
  for k = 1:numel(use), if i ~= k && ~isempty(strfind(char(X{k}.' + 48), a)), hits = hits + 1; fprintf('  %s inside %s\n', use(i), use(k)); end, end, end
fprintf('cross-file first-day matches %d\n', hits);
% phase check: mean activity by 30-min bin of the 24 h day, all flies, all whole days
prof = zeros(48,1);
for i = 1:numel(use), x = X{i}; nd = floor(numel(x)/86400); x = reshape(x(1:nd*86400), 1800, 48, nd); prof = prof + squeeze(mean(sum(x,1),3)).' / numel(use); end
[~, o] = sort(prof, 'descend');
fprintf('top 6 half-hour bins (start hour): %s\n', mat2str((o(1:6)-1)/2));
fprintf('mean counts/30 min  ZT0-12: %.1f   ZT12-24: %.1f\n', mean(prof(1:24)), mean(prof(25:48)));
fprintf('bins around ZT0: %s | around ZT12: %s\n', mat2str(round(prof([47 48 1 2 3])')), mat2str(round(prof(23:27)')));
tbl = table(use(:), n, n/86400, fp, 'VariableNames', {'File','Samples','Days','SHA256prefix'});
writetable(tbl, fullfile(fileparts(T), 'csld_manifest.csv'));
for i = 1:height(tbl), fprintf('%-50s %5.2f %s\n', tbl.File(i), tbl.Days(i), tbl.SHA256prefix(i)); end
