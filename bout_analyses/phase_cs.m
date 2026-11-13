T = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\csld_stage';
f = dir(fullfile(T, '*.txt')); names = string({f.name}); exps = extractBefore(names, regexpPattern('_(multibeam|mb|4b|2boardmb)'));
for e = unique(exps)
  idx = find(exps == e); prof = zeros(96,1);
  for i = idx, x = readmatrix(fullfile(T, names(i)), 'FileType', 'text'); x = x(:,1); x = x(~isnan(x));
    nd = floor(numel(x)/86400); x = reshape(x(1:nd*86400), 900, 96, nd); prof = prof + squeeze(mean(sum(x,1),3)).'; end
  [~, o] = max(prof(41:56)); [~, o0] = max(prof([89:96 1:8]));
  fprintf('%-10s flies %2d  largest 15-min bin near ZT12 starts at ZT %.2f; near ZT0 at ZT %.2f\n', e, numel(idx), (40+o-1)/4, mod((88+o0-1)/4, 24));
end
