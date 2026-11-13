base = 'C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data';
scr  = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad';
% CS DD
T = fullfile(scr, 'csdd_stage'); if ~isfolder(T), mkdir(T); end
f = dir(fullfile(base, 'CS DD', '*_from1stLD.txt'));
for i = 1:numel(f), copyfile(fullfile(f(i).folder, f(i).name), fullfile(T, f(i).name)); end
fprintf('CS DD staged: %d\n', numel(dir(fullfile(T, '*.txt'))));
% per0 LD, one folder per experiment
ex = struct('name', {'8_18_2017','9_21_2020','8_25_2017','9_22_2017','9_29_2017'}, ...
            'startZT', {11.0, 11.5, 20.0, 20.5, 8.0}, 'Telim', {13, 12.5, 4, 3.5, 16});
f = dir(fullfile(base, 'per0 ld', '*.txt')); names = string({f.name});
for e = ex
    T = fullfile(scr, 'per0ld_stage', e.name); if ~isfolder(T), mkdir(T); end
    m = names(startsWith(names, e.name + "_"));
    for k = m, copyfile(fullfile(base, 'per0 ld', k), fullfile(T, k)); end
    % day-by-day transition check on the trimmed record (ZT0 at sample 1)
    P = [];
    for k = m
        x = readmatrix(fullfile(T, k), 'FileType', 'text'); x = x(:, 1); x = x(~isnan(x));
        x = x(round(e.Telim * 3600) + 1:end);
        nd = floor(numel(x) / 86400);
        P(:, :, end+1) = squeeze(sum(reshape(x(1:nd*86400), 900, 96, nd), 1)); %#ok<SAGROW>
    end
    P = mean(P(:, :, 2:end), 3);               % 96 x nd, fly mean
    nd = size(P, 2); ok = false(nd, 2);
    for d = 1:nd
        y = P(:, d); prev = [P(end, max(d-1, 1)); y(1:end-1)];
        rise = y - prev;                        % rise INTO each 15-min bin
        thr = prctile(rise, 90);
        w0 = [93:96 1:5]; w12 = 45:53;          % +-1 h around ZT0 and ZT12 (bins start at ZT (b-1)/4)
        if d == 1, w0 = 1:5; end
        ok(d, :) = [max(rise(w0)) >= thr, max(rise(w12)) >= thr];
    end
    fprintf('%-10s flies %2d, days %d: lights-on rise found on %d/%d days, lights-off on %d/%d\n', ...
        e.name, numel(m), nd, sum(ok(:, 1)), nd, sum(ok(:, 2)), nd);
end
