base = 'C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data';
sets = {'CS DD','per0 ld','per0 dd','CS ld'};
md = java.security.MessageDigest.getInstance('SHA-256');
allHash = strings(0); allName = strings(0); firstDay = {};
for s = 1:numel(sets)
  f = dir(fullfile(base, sets{s}, '*.txt')); names = string({f.name});
  if s == 3   % per0 dd: the 48 used
    cut = names(endsWith(names, '_from1stLD.txt')); names = names(~ismember(names, erase(cut, '_from1stLD.txt')));
  end
  fprintf('\n== %s: %d files\n', sets{s}, numel(names));
  prof = zeros(96, 1); days = zeros(numel(names), 1); hp = strings(numel(names), 1);
  for i = 1:numel(names)
    p = fullfile(base, sets{s}, names(i));
    fid = fopen(p); b = fread(fid, inf, '*uint8'); fclose(fid);
    md.reset(); h = typecast(md.digest(b), 'uint8'); hp(i) = string(reshape(dec2hex(h(1:8))', 1, []));
    x = readmatrix(p, 'FileType', 'text'); x = x(:, 1); x = x(~isnan(x));
    days(i) = numel(x) / 86400;
    allName(end+1) = sets{s} + "/" + names(i); firstDay{end+1} = char(x(1:min(end, 86400)).' + 48); %#ok<SAGROW>
    md.reset(); h = typecast(md.digest(typecast(int32(x), 'uint8')), 'uint8'); allHash(end+1) = string(reshape(dec2hex(h(1:8))', 1, [])); %#ok<SAGROW>
    if s <= 2
      nd = floor(numel(x) / 86400); y = reshape(x(1:nd*86400), 900, 96, nd); prof = prof + squeeze(mean(sum(y, 1), 3)).';
      if s == 2
        q = reshape(x(1:nd*86400), 900, 96, nd); pq = squeeze(mean(sum(q, 1), 3)); [~, o] = sort(pq, 'descend');
        fprintf('  %-38s %5.2f d  top 15-min bins start at h %s\n', names(i), days(i), mat2str((o(1:3) - 1) / 4));
      end
    end
  end
  if s ~= 2, fprintf('  days: %s\n', mat2str(round(days.', 2))); end
  if s <= 2
    [~, o] = sort(prof, 'descend'); fprintf('  pooled: top 15-min bins start at hours %s\n', mat2str((o(1:6) - 1).' / 4));
  end
end
[u, ~, j] = unique(allHash); c = accumarray(j, 1);
fprintf('\ndistinct recordings across all four sets: %d of %d\n', numel(u), numel(allHash));
hits = 0;
for i = 1:numel(firstDay), for k = 1:numel(firstDay)
  if i ~= k && numel(firstDay{i}) == 86400 && contains(firstDay{k}, firstDay{i}), hits = hits + 1; fprintf('  first day of %s inside %s\n', allName(i), allName(k)); end
end, end
fprintf('cross-set first-day matches: %d\n', hits);
