function test_BeamFingerprint(rawFile)
%TEST_BEAMFINGERPRINT  Golden test: BeamFingerprint's preprocessing must reproduce beam_hmm/pipeline.py
%   on week 1 (golden_week1_qc.json, written by the Python pipeline). Usage, from beam_hmm/matlab:
%       test_BeamFingerprint('path/to/raw_positions_all64_week1.txt')
%   Checks, channels 17-32: spikes fixed, hub value/share/returns, samples replaced, position changes
%   (total and first 24 h), share of >= 4-beam jumps after repair. Integers must match exactly; shares to 1e-12.
assert(strcmp(fileparts(which('BeamFingerprint')), pwd), 'run from beam_hmm/matlab (a stale copy would win otherwise)');
G = jsondecode(fileread('golden_week1_qc.json'));
out = BeamFingerprint(rawFile, 'Channels', 17:32, 'Plot', false);
Q = out.QC; nfail = 0;
assert(numel(out.Time) == G.n_samples, 'sample count %d vs golden %d', numel(out.Time), G.n_samples);
for i = 1:numel(G.channels)
    if iscell(G.channels), g = G.channels{i}; else, g = G.channels(i); end
    k = find(Q.Channel == g.channel);
    chk = {'SpikesFixed', g.spikes; 'HubSamplesReplaced', g.replaced; 'Changes', g.changes; 'Changes24h', g.changes24};
    for r = 1:size(chk, 1)
        if Q.(chk{r, 1})(k) ~= chk{r, 2}
            fprintf('[FAIL] ch%d %s: %g vs golden %g\n', g.channel, chk{r, 1}, Q.(chk{r, 1})(k), chk{r, 2}); nfail = nfail + 1;
        end
    end
    if ~isempty(g.hub) && Q.HubValue(k) ~= g.hub
        fprintf('[FAIL] ch%d hub %g vs %g\n', g.channel, Q.HubValue(k), g.hub); nfail = nfail + 1;
    end
    fl = {'HubShare', g.share; 'HubReturns', g.returns; 'FracJumpsGE4_final', g.frac_big_final};
    for r = 1:size(fl, 1)
        if ~isempty(fl{r, 2}) && abs(Q.(fl{r, 1})(k) - fl{r, 2}) > 1e-12
            fprintf('[FAIL] ch%d %s: %.15g vs %.15g\n', g.channel, fl{r, 1}, Q.(fl{r, 1})(k), fl{r, 2}); nfail = nfail + 1;
        end
    end
end
if nfail == 0
    fprintf('test_BeamFingerprint: all %d channels match the Python pipeline\n', numel(G.channels));
else
    error('test_BeamFingerprint: %d mismatches', nfail);
end
end
