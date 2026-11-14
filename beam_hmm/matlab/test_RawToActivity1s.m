function test_RawToActivity1s(rawFile, goldenFile)
%TEST_RAWTOACTIVITY1S  Golden test of RawToActivity1s against the Python reference (beam_hmm/labconv.py).
%   test_RawToActivity1s('path/to/raw_positions_all64_week1.txt')        % uses golden_week1_activity.json
%   Run from beam_hmm/matlab. For channels 17-32 it checks four variants -- legacy (PreProcess +
%   multibeam_locomotion), beam_hmm (repairs, 0 kept), beam_hmm with 0 folded, and beam_hmm without repairs --
%   on four numbers each: length, total beams, seconds with movement, and a position-weighted checksum
%   sum(v .* mod(1:n, 9973)). All must match exactly.
if nargin < 2, goldenFile = 'golden_week1_activity.json'; end
assert(strcmp(fileparts(which('RawToActivity1s')), pwd), 'run from beam_hmm/matlab (a stale copy would win otherwise)');
G = jsondecode(fileread(goldenFile));
runs = {'legacy', {'Mode', 'legacy'}; 'beam_hmm', {}; 'beam_hmm_fold0', {'KeepBeyondBeam1', false}; 'norepair', {'Repair', false}};
nfail = 0;
for r = 1:size(runs, 1)
    o = RawToActivity1s(rawFile, 'Channels', 17:32, runs{r, 2}{:});
    for k = 1:16
        c = 16 + k; g = G.channels.(sprintf('x%d', c)).(runs{r, 1});
        v = o.Counts(:, k); n = numel(v); w = mod((1:n)', 9973);
        got = [n, sum(v), sum(v > 0), sum(v .* w)]; want = [g.n, g.sum, g.nonzero, g.wsum];
        if any(got ~= want)
            fprintf('[FAIL] %s ch%d: got [%d %d %d %d], golden [%d %d %d %d]\n', runs{r, 1}, c, got, want); nfail = nfail + 1;
        end
    end
end
if nfail == 0
    fprintf('test_RawToActivity1s: all 4 variants x 16 channels match the Python reference\n');
else
    error('test_RawToActivity1s: %d mismatches', nfail);
end
end
