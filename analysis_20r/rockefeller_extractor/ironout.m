function [activity_out, nBridged, nBridgedSamples] = ironout(activity_wave, maxWakeSeconds, timeInterval)
%IRONOUT Absorb brief wake interruptions that sit between two sleep bouts.
%
%   activity_out = IRONOUT(activity_wave, maxWakeSeconds, timeInterval)
%   relabels as sleep (0) any wake run (1s) lasting no more than
%   maxWakeSeconds that has sleep on BOTH sides, merging the two flanking
%   sleep bouts into one. The biological argument is that an isolated beam
%   break inside a long quiescent stretch is more likely a postural adjustment
%   than a genuine arousal.
%
%   [activity_out, nBridged, nBridgedSamples] = IRONOUT(...) also returns how
%   many wake runs were absorbed and how many samples that amounted to, so the
%   size of the effect can be measured rather than assumed.
%
%   CONVENTION
%   ----------
%   activity_wave is a DOUBLE vector with 0 = SLEEP and 1 = WAKE -- the
%   activity wave returned by sleeptrace, and the convention used throughout
%   this pipeline. The output uses the same convention, so sleeptrace output
%   can be passed straight in:
%
%       wave = sleeptrace(activity, dt);
%       wave = ironout(wave, 1, dt);
%
%   A LOGICAL input is refused. An earlier version of this codebase used
%   logical true = SLEEP, the exact opposite. A stale mask of that kind would
%   otherwise be accepted and bridged backwards -- brief SLEEP runs absorbed
%   into wake -- silently and without error. If your data are genuinely
%   0 = sleep, 1 = wake but happen to be logical, wrap them in double().
%
%   WHAT THIS CHANGES
%   -----------------
%   This is a deliberate departure from the conventional five-minute inactivity
%   rule. Samples in which the animal actually moved are scored as sleep.
%
%     * Total sleep rises slightly; total wake falls by the same amount.
%     * Sleep bouts get longer and fewer, because separate bouts are merged.
%     * THE WAKE-BOUT DISTRIBUTION IS TRUNCATED AT ITS SHORT END. Every
%       interior wake run at or below the threshold is removed by construction,
%       so the shortest wake bouts no longer exist downstream. Anything fitted
%       to wake-bout durations - exponential, power-law, hyperexponential, HMM
%       dwell times - is affected, and such fits are NOT comparable to fits
%       from unbridged data. When modelling bout-length distributions, prefer
%       maxWakeSeconds = 0 and handle brief arousals in the model instead.
%
%   The operation cannot create a sleep bout: both flanks are already sleep, so
%   an existing bout can only be extended. Every sleep bout therefore still
%   satisfies the five-minute rule afterwards.
%
%   BOUNDARY CONVENTION
%   -------------------
%   Only INTERIOR wake runs are bridged. A run touching the first or last
%   sample has sleep on one side only, so there is nothing to merge, and its
%   true duration is censored by the recording boundary in any case.
%
%   All qualifying runs are filled from the ORIGINAL wave in a single pass, so
%   the result does not depend on evaluation order. The operation is
%   idempotent: applying it twice with the same threshold changes nothing the
%   second time, because bridging only lengthens sleep runs and never creates a
%   new short wake run.
%
%   INPUTS
%   ------
%   activity_wave   Nonempty double column vector of 0s and 1s,
%                   0 = sleep, 1 = wake.
%   maxWakeSeconds  Longest wake run to absorb, in seconds. Nonnegative. Zero
%                   disables bridging and returns the input unchanged.
%   timeInterval    Duration of one sample, in seconds. Positive.
%
%   The threshold is converted with floor(maxWakeSeconds / timeInterval), so a
%   run of N samples qualifies when N*timeInterval <= maxWakeSeconds. Coarse
%   sampling therefore disables bridging automatically: at timeInterval = 60 s
%   no run can be 1 s long, and the input is returned unchanged.
%
%   EXAMPLE
%   -------
%       % 1 s bins: a lone active second inside a long quiet stretch
%       w = [zeros(600,1); 1; zeros(600,1)];
%       [w2, n] = ironout(w, 1, 1);       % n = 1, w2 is all zeros
%
%       % the same gap is left alone when it is two seconds long
%       w = [zeros(600,1); ones(2,1); zeros(600,1)];
%       [~, n] = ironout(w, 1, 1);        % n = 0
%
%   See also SLEEPTRACE, SEGMENT_SLEEP_PHASES, PROCESS_LOCOMTOSLEEP

arguments
    % No class declared on purpose: declaring 'double' would make MATLAB
    % convert a logical input silently BEFORE validation, defeating the check.
    activity_wave (:,1) {mustBeNonempty, mustBeActivityWave}
    maxWakeSeconds (1,1) {mustBeNumeric, mustBeReal, mustBeNonnegative, mustBeFinite} = 1
    timeInterval (1,1) {mustBeNumeric, mustBeReal, mustBePositive, mustBeFinite} = 1
end

activity_out = activity_wave;
nBridged = 0;
nBridgedSamples = 0;

maxBins = floor(maxWakeSeconds / timeInterval);
if maxBins < 1
    return
end

% Locate every contiguous run of wake (1s).
isWake = (activity_out == 1);
transitions = diff([false; isWake; false]);
starts = find(transitions == 1);
ends   = find(transitions == -1) - 1;
lengths = ends - starts + 1;

% Interior runs only, and only those at or below the threshold.
interior = starts > 1 & ends < numel(activity_out);
qualify  = interior & lengths <= maxBins;
if ~any(qualify)
    return
end

for k = find(qualify).'
    activity_out(starts(k):ends(k)) = 0;
end

nBridged = nnz(qualify);
nBridgedSamples = sum(lengths(qualify));

end

%% ------------------------------------------------------------------------
function mustBeActivityWave(x)
%MUSTBEACTIVITYWAVE Double vector of 0s and 1s, 0 = sleep, 1 = wake.
if ~isa(x, 'double')
    error('ironout:NotActivityWave', ...
        ['Input must be a DOUBLE activity wave (0 = sleep, 1 = wake), as ' ...
         'returned by sleeptrace; got %s. Logical input is refused because ' ...
         'the earlier convention here was logical true = SLEEP, and a stale ' ...
         'mask would be bridged backwards without error. If your data are ' ...
         'already 0 = sleep, 1 = wake, pass double(x).'], class(x));
end
if any(x ~= 0 & x ~= 1)          % NaN also fails this test
    error('ironout:NotBinary', ...
        'Activity wave must contain only 0 (sleep) and 1 (wake).');
end
end