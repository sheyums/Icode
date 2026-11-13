[X, ds] = qload(); fprintf('flies %d; lag pairs fly1: %d %d %d\n', numel(X), numel(cat(1, X{1}.P{1}{:})), numel(cat(1, X{1}.P{2}{:})), numel(cat(1, X{1}.P{3}{:})));
tic; rs = RandStream('mrg32k3a', 'Seed', 1);
for sc = 1:3
    Y = cellfun(@(x) qsim(x, sc, rs), X, 'UniformOutput', false);
    v = cell2mat(cellfun(@(y) qstats(y.d, y), Y, 'UniformOutput', false).');
    fprintf('sc %d raw mean r1 %.3f r2 %.3f r3 %.3f\n', sc, mean(v(:, 1:3), 'omitnan'));
end
R = qrun(X(1:10), 20, 1); disp(R.F); [k, ci] = qkappa(R.fast, ds(1:10), 50, 1); fprintf('kappa %.3f [%.3f %.3f]; 10 flies x 20 perms %.1f s\n', k, ci, toc);
