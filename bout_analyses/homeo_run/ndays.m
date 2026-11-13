L = load('homeo_results.mat', 'flies');
for s = ["per0DD", "CSLD"]
    F = L.flies.(s); d = arrayfun(@(f) (f.Onset(end) + f.Dur(end) - f.Onset(1)) / 86400, F); ns = arrayfun(@(f) sum(f.IsSleep), F);
    fprintf('%s: %d flies; days min %.2f median %.2f max %.2f; whole days 0..9 %s; sleep bouts/fly median %g\n', s, numel(F), min(d), median(d), max(d), mat2str(histcounts(floor(d), 0:10)), median(ns));
end
