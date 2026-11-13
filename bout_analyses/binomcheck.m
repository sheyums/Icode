src = fileread('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode\BoutIrreversibility.m');
k = strfind(src, 'function p = binom_upper'); body = src(k:end);
d = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\bc'; if ~isfolder(d), mkdir(d); end
fid = fopen(fullfile(d, 'binom_upper.m'), 'w'); fwrite(fid, body); fclose(fid); addpath(d);
mx = 0;
for n = [1 3 10 40 48 200], for kk = 0:n
  mx = max(mx, abs(binom_upper(kk, n, 0.05) - (1 - binocdf(kk - 1, n, 0.05))));
end, end
fprintf('max abs difference from 1 - binocdf: %.3g\n', mx);
