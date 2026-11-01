function sleep_out = sleeptrace(activity,timeInterval)
%SLEEPTRACE Convert locomotor activity to a binary sleep-wake time series.
%
%   sleep_out = SLEEPTRACE(activity,timeInterval) classifies each sample of
%   a locomotor activity time series as sleep or wake using the conventional
%   retrospective five-minute inactivity rule:
%
%       sleep_out = 0  -> sleep
%       sleep_out = 1  -> wake
%
%   An inactivity bout is defined as a contiguous run of samples for which
%   activity <= cutoff. If the complete inactivity bout lasts at least five
%   minutes, every sample in that bout is classified as sleep. Thus, this
%   function uses retrospective whole-bout scoring, also called
%   Interpretation B.
%
%   INPUTS
%   ------
%   activity
%       Nonempty real numeric vector containing locomotor activity values.
%       The vector may be a row or column vector. Values at or below cutoff
%       are treated as inactive.
%
%       Values must be finite; NaN and Inf are rejected by the arguments
%       block and raise an error. Negative values are also invalid, but are
%       deliberately NOT rejected by a validator, so that a batch caller is
%       not terminated by one bad record. See INVALID ACTIVITY VALUES below.
%
%   timeInterval
%       Positive finite scalar specifying the duration represented by each
%       activity sample, in seconds. For example:
%
%           timeInterval = 1;    % one-second activity bins
%           timeInterval = 60;   % one-minute activity bins
%
%   OUTPUT
%   ------
%   sleep_out
%       Double-precision binary vector having the same size and orientation
%       as activity. Zero represents sleep and one represents wake.
%
%       This holds only when the trace was classified. If it was refused,
%       sleep_out is a scalar NaN irrespective of the size of activity. See
%       INVALID ACTIVITY VALUES below.
%
%   INVALID ACTIVITY VALUES
%   -----------------------
%   If activity contains any negative value the trace is NOT classified.
%   A warning with identifier
%
%       sleeptrace:NegativeActivity
%
%   is issued and sleep_out is returned as a scalar NaN. This is a refusal,
%   not a partial result: no sample is classified, even where an inactivity
%   run elsewhere in the trace would otherwise have qualified as sleep.
%
%   The refusal is a warning rather than an error so that a caller looping
%   over many files can skip one bad record and continue. Callers must
%   therefore check the return value before using it, because a scalar NaN
%   will otherwise propagate silently, and 0/1 doubles cast to logical
%   without complaint:
%
%       sleep_out = sleeptrace(activity,timeInterval);
%       if isscalar(sleep_out) && isnan(sleep_out)
%           % skip this record
%       end
%
%   NOTE ON POLARITY
%   ----------------
%   sleep_out is an ACTIVITY wave: one marks wake. Downstream code that
%   expects a SLEEP mask, notably segment_sleep_phases, requires the
%   opposite convention and must invert first:
%
%       sleep_mask = (sleeptrace(activity,timeInterval) == 0);
%
%   Because a vector of 0s and 1s is a perfectly legal logical, MATLAB will
%   cast it without any warning, and omitting the inversion yields the exact
%   complement of every downstream statistic rather than an error.
%
%   DEFINITION AND BOUNDARY CONVENTION
%   ----------------------------------
%   The inactivity cutoff and minimum sleep duration are defined internally:
%
%       cutoff       = 0;
%       sleepDuration = 5 minutes;
%
%   The duration of an inactivity bout containing N samples is calculated as
%
%       boutDuration = N*timeInterval.
%
%   Therefore, a bout qualifies as sleep when
%
%       N*timeInterval >= 300 seconds.
%
%   Once a bout qualifies, the entire inactivity bout, including its first
%   five minutes, is labeled as sleep.
%
%   EXAMPLE
%   -------
%   For one-minute activity bins:
%
%       activity     = [3 0 0 0 0 0 0 2];
%       sleep_out    = sleeptrace(activity,60);
%
%   The result is:
%
%       sleep_out =
%            1     0     0     0     0     0     0     1
%
%   The six-minute inactivity bout is classified entirely as sleep.
%
%   See also SEGMENT_SLEEP_PHASES, PROCESS_LOCOMTOSLEEP, STRESSTEST_SLEEPTRACE
%
%% Took an earlier version and made it look smarter with help from GPT 5.6. Sheyum August 6, 2026
%%

arguments
    activity {mustBeNumeric,mustBeReal,mustBeVector,mustBeNonempty,...
        mustBeFinite}
    timeInterval (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
end

% Negative locomotor activity is invalid. Return a scalar NaN rather than
% terminating the caller with an error.
if any(activity < 0,"all")
    numberNegative = nnz(activity < 0);
    firstNegative = find(activity < 0,1,"first");

    % NOTE: the message must be built with single-quoted char concatenation.
    % Bracketing double-quoted literals, ["a" "b"], builds a 1-by-2 string
    % ARRAY rather than one string, which warning() rejects outright.
    warning('sleeptrace:NegativeActivity', ...
        ['Activity contains %d negative value(s). The first occurs at ' ...
        'sample %d with value %.15g. Sleep classification was not ' ...
        'performed, and sleep_out is returned as NaN.'], ...
        numberNegative,firstNegative,activity(firstNegative));

    sleep_out = NaN;
    return
end

% Internal classification parameters.
cutoff = 0;
sleepDuration = 5*60;   % Minimum sleep-bout duration, in seconds.

% Initially classify every sample as wake.
sleep_out = ones(size(activity));

% Identify samples satisfying the inactivity criterion.
inactive = activity <= cutoff;

% Add false boundary values so inactivity bouts beginning at the first
% sample or ending at the final sample are handled identically to internal
% bouts.
transitions = diff([false; inactive(:); false]);

% Locate the first and last sample of every contiguous inactivity bout.
boutStarts = find(transitions == 1);
boutEnds   = find(transitions == -1) - 1;

% Determine which entire inactivity bouts satisfy the five-minute rule.
boutDurations = (boutEnds - boutStarts + 1)*timeInterval;
isSleepBout = boutDurations >= sleepDuration;

% Retrospectively classify every sample in each qualifying bout as sleep.
for k = find(isSleepBout).'
    sleep_out(boutStarts(k):boutEnds(k)) = 0;
end

end
