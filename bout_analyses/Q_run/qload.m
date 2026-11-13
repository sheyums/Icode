function [X, ds] = qload()
L = load(fullfile('..', 'homeo_run', 'homeo_results.mat'), 'flies');
sets = {'per0DD', 24; 'CSLD', 0; 'CSDD', 24}; X = {}; ds = [];
for k = 1:3
    F = L.flies.(sets{k, 1});
    for f = 1:numel(F), X{end+1} = qprep(F(f), sets{k, 2}); ds(end+1) = k; end %#ok<AGROW>
end
end
