function segmented_out = segment_sleep_phases(activity_wave, options)
% SEGMENT_SLEEP_PHASES Separates a binary activity wave into phase components and calculates bout architecture.
%
%   This function processes a binary sleep-wake time series to facilitate
%   two distinct types of downstream circadian analysis:
%
%   1. Total Sleep Time (Hard-Masking): Uses a strict phase mask to extract
%      exact phase bins (e.g., Light vs. Dark). This ensures every minute
%      of sleep is credited to the correct phase, even if a bout crosses a
%      phase boundary.
%   2. Sleep Architecture (Bout Extraction): Identifies the start, end, and
%      duration of every intact bout. It assigns an "Origin Phase" based on
%      when the bout began, preventing artificial bout-fracturing at phase
%      boundaries. Boundary-censored bouts (touching the start/end of the
%      recording) are automatically excluded from architecture stats.
%
%   Syntax:
%       segmented_out = segment_sleep_phases(activity_wave)
%       segmented_out = segment_sleep_phases(activity_wave, Name=Value, ...)
%
%   Input Arguments:
%       activity_wave - (double) N-by-1 vector of 0s and 1s,
%                       0 = SLEEP, 1 = WAKE. This is exactly what sleeptrace
%                       returns, so its output can be passed straight in:
%
%                           wave = sleeptrace(activity, timeInterval);
%                           out  = segment_sleep_phases(wave, ...);
%
%                       A LOGICAL input is refused. An earlier version of this
%                       function took logical true = SLEEP, the opposite
%                       convention; accepting logical would let a stale mask of
%                       that kind through with .sleep_bouts and .wake_bouts
%                       swapped and every total complemented, silently. If your
%                       data are genuinely 0 = sleep but logical, pass double().
%
%   Name-Value Arguments:
%       mode         - (numeric) Segmentation mode.
%                      1: Entire trace only (No phase splitting). (Default)
%                      2: Light / Dark split.
%                      3: Subjective Day / Subjective Night (Constant Darkness).
%       timeInterval - (numeric) Duration of each sample bin in seconds. (Default: 1)
%       dayLength    - (numeric) Total hours in one complete circadian cycle. (Default: 24)
%       daytime      - (numeric) Duration of the light (or subjective day) phase in hours.
%                      Must not exceed 'dayLength'. (Default: 12)
%       nighttime    - (numeric) Duration of the dark (or subjective night) phase in
%                      hours. Optional consistency check only: the dark phase is always
%                      derived as dayLength - daytime, so supplying 'nighttime' merely
%                      asserts that daytime + nighttime == dayLength and errors if it
%                      does not. (Default: [], meaning no check)
%       phaseOffset  - (numeric) The ZT/CT hour at which the trimmed activity
%                      wave begins. Used to align the phase mask. (Default: 0)
%       verbose      - (logical) Print a message when boundary-censored bouts are
%                      dropped. The counts are always available in the output as
%                      .stats.*.n_censored regardless of this flag. (Default: false)
%
%   Output Structure (segmented_out):
%       All binary arrays use the input convention: DOUBLE, 0 = sleep, 1 = wake.
%
%       .entire       - The original, unbroken activity wave.
%       .light        - (Mode 2) All Light phase bins.
%       .dark         - (Mode 2) All Dark phase bins.
%       .subj_day     - (Mode 3) All Subjective Day bins.
%       .subj_night   - (Mode 3) All Subjective Night bins.
%       .stats        - Structure containing architecture data for:
%           .sleep_bouts    (runs of 0s)
%           .wake_bouts     (runs of 1s)
%               .starts       - Array index where each bout begins.
%               .ends         - Array index where each bout ends.
%               .durations    - Length of each bout in *DATA POINTS*.
%                               (Multiply by timeInterval to get seconds).
%               .origin_phase - String array indicating which phase the bout started in.
%               .n_censored   - Number of boundary-censored bouts that were excluded
%                               (0, 1, or 2). Useful for sanity-checking short traces.
%
%   NOTE ON .light / .dark: these are concatenations of non-contiguous bins, so
%   neighbouring elements are not neighbours in time. They are correct for TOTAL
%   sleep time per phase, but running any bout/run-length analysis on them will
%   fabricate bouts at every phase boundary. Use .stats with .origin_phase for
%   per-phase bout architecture instead.
%
%   Do NOT assume numel(.light) == numel(.dark). They are equal only when the
%   trace spans a whole number of cycles AND the photoperiod is symmetric. A
%   partial final cycle, or a nonzero phaseOffset, makes them differ. For rates
%   and fractions divide by each phase's own bin count, which mean() does:
%
%       mean(out.light == 0)   % sleep FRACTION during light, normalized correctly
%       mean(out.light)        % WAKE fraction during light (1 = wake)
%       sum(out.light == 0)    % sleep BINS during light, a total -- do not normalize
%
%   Reconstructing a denominator by hand (e.g. assuming a 50/50 split, or
%   dayLength*bins_per_hour/2) is the failure mode this warns against.
%
%   Example:
%       % Segment a 24h LD cycle starting at ZT0, with 60-second bins
%       wave = sleeptrace(activity, 60);
%       sleep_data = segment_sleep_phases(wave, mode=2, timeInterval=60);
%
%       % Total daytime sleep in minutes. sum(... == 0) counts SLEEP BINS, so
%       % scale by the bin duration: minutes = bins * timeInterval / 60.
%       day_sleep_min = sum(sleep_data.light == 0) * (60 / 60);   % 60 s bins -> x1
%
%       % Find the duration (in data points) of the first sleep bout
%       first_bout_length = sleep_data.stats.sleep_bouts.durations(1);
%
%       % Per-phase bout architecture. Select bouts by ORIGIN PHASE; do not run
%       % bout analysis on .light / .dark (see the note above).
%       b = sleep_data.stats.sleep_bouts;
%       light_bout_durations = b.durations(b.origin_phase == "Light");
%
%   Errors:
%       segment_sleep_phases:NotActivityWave       input not double
%       segment_sleep_phases:NotBinary             input not all 0/1
%       segment_sleep_phases:DaytimeExceedsCycle   daytime > dayLength
%       segment_sleep_phases:PhaseDurationMismatch daytime + nighttime ~= dayLength
%       segment_sleep_phases:NighttimeNotScalar    nighttime is non-scalar
%
%   See also SLEEPTRACE, PROCESS_LOCOMTOSLEEP, CALCBOUTSIZES

arguments
    % No class declared on purpose: declaring 'double' would make MATLAB
    % convert a logical input silently BEFORE validation, defeating the check.
    activity_wave (:,1) {mustBeActivityWave}
    options.mode (1,1) {mustBeMember(options.mode, [1, 2, 3])} = 1
    options.timeInterval (1,1) {mustBeNumeric, mustBePositive} = 1
    options.dayLength (1,1) {mustBeNumeric, mustBePositive} = 24
    options.daytime (1,1) {mustBeNumeric, mustBeNonnegative} = 12
    options.nighttime {mustBeNumeric, mustBeNonnegative} = []
    options.phaseOffset (1,1) {mustBeNumeric, mustBeNonnegative} = 0
    options.verbose (1,1) logical = false
end

SLEEP = 0;
WAKE  = 1;

% Calculate bin metrics for the full cycle
bins_per_hour    = 3600 / options.timeInterval;
total_cycle_bins = round(options.dayLength * bins_per_hour);
day_bins         = round(options.daytime * bins_per_hour);

% A daytime longer than the cycle would make night_bins negative. false(-n,1)
% silently yields an empty array, leaving cycle_template shorter than one full
% cycle, which corrupts the period of the replicated mask without any error.
% Reject it up front instead.
if day_bins > total_cycle_bins
    error('segment_sleep_phases:DaytimeExceedsCycle', ...
        'daytime (%g hrs) cannot exceed dayLength (%g hrs).', ...
        options.daytime, options.dayLength);
end

% The dark phase is always the remainder of the cycle. If the caller also
% supplied nighttime, assert the three durations are mutually consistent
% rather than silently ignoring it.
night_bins = total_cycle_bins - day_bins;
if ~isempty(options.nighttime)
    if ~isscalar(options.nighttime)
        error('segment_sleep_phases:NighttimeNotScalar', ...
            'nighttime must be a scalar or empty.');
    end
    if round(options.nighttime * bins_per_hour) ~= night_bins
        error('segment_sleep_phases:PhaseDurationMismatch', ...
            ['daytime (%g hrs) + nighttime (%g hrs) = %g hrs, which does not ' ...
            'match dayLength (%g hrs).'], ...
            options.daytime, options.nighttime, ...
            options.daytime + options.nighttime, options.dayLength);
    end
end

% Base output always contains the unbroken trace
segmented_out = struct();
segmented_out.entire = activity_wave;

% Build PHASE mask (true = Light/SubjDay, false = Dark/SubjNight). This is a
% phase indicator, unrelated to the sleep/wake coding of the data.
cycle_template = [true(day_bins, 1); false(night_bins, 1)];

% Shift mask if an offset exists
offset_bins = round(options.phaseOffset * bins_per_hour);
if offset_bins > 0
    shift_amount = mod(offset_bins, total_cycle_bins);
    cycle_template = circshift(cycle_template, -shift_amount);
end

% Replicate mask to match the length of the activity wave
num_cycles = ceil(length(activity_wave) / total_cycle_bins);
full_mask  = repmat(cycle_template, num_cycles, 1);
full_mask  = full_mask(1:length(activity_wave));

% 1. Hard-Boundary Segmentation (Perfect for calculating Total Sleep Time per phase)
if options.mode == 2
    % Mode 2: Light / Dark Split
    segmented_out.light = activity_wave(full_mask);
    segmented_out.dark  = activity_wave(~full_mask);
elseif options.mode == 3
    % Mode 3: Subjective Day / Subjective Night Split
    segmented_out.subj_day   = activity_wave(full_mask);
    segmented_out.subj_night = activity_wave(~full_mask);
end

% 2. Sleep/Wake Architecture Analysis (Perfect for Bout number & length)
[s_starts, s_ends, s_durs, s_origins, s_ncens] = ...
    extract_bouts(activity_wave, SLEEP, full_mask, options.mode, options.verbose);
[w_starts, w_ends, w_durs, w_origins, w_ncens] = ...
    extract_bouts(activity_wave, WAKE,  full_mask, options.mode, options.verbose);

segmented_out.stats.sleep_bouts.starts       = s_starts;
segmented_out.stats.sleep_bouts.ends         = s_ends;
segmented_out.stats.sleep_bouts.durations    = s_durs;
segmented_out.stats.sleep_bouts.origin_phase = s_origins;
segmented_out.stats.sleep_bouts.n_censored   = s_ncens;

segmented_out.stats.wake_bouts.starts        = w_starts;
segmented_out.stats.wake_bouts.ends          = w_ends;
segmented_out.stats.wake_bouts.durations     = w_durs;
segmented_out.stats.wake_bouts.origin_phase  = w_origins;
segmented_out.stats.wake_bouts.n_censored    = w_ncens;

end

%% Local Helper Functions
function [starts, ends, durations, origins, nCensored] = ...
    extract_bouts(timeseries, target_state, mask, mode, verbose)
% Extracts every run of 'target_state' (0 = sleep, 1 = wake) and applies
% standard boundary-censoring logic.

x = double(timeseries);
x(timeseries ~= target_state) = 100;
x(timeseries == target_state) = 101;
x = [100; x(:); 100];

% Find bout starts and calculate consecutive sequence lengths
uptodownlocs = find(diff(x) == 1) + 1;
length_s = find(diff(x) == -1) + 1 - uptodownlocs;

% Derive start and end indices relative to the original unpadded timeseries
starts = uptodownlocs - 1;
ends = starts + length_s - 1;

% Identify and exclude boundary-censored bouts touching the start or end
isCensored = false(size(length_s));
if ~isempty(length_s)
    if timeseries(1) == target_state
        isCensored(1) = true;
    end
    if timeseries(end) == target_state
        isCensored(end) = true;
    end
end

nCensored = nnz(isCensored);
if verbose && nCensored > 0
    state_name = 'Sleep';
    if target_state == 1
        state_name = 'Wake';
    end
    fprintf('segment_sleep_phases (%s): %d of %d bout(s) touch the start/end of the trace and were excluded as censored.\n', ...
        state_name, nCensored, numel(length_s));
end

% Remove censored bouts from all tracked arrays
starts(isCensored) = [];
ends(isCensored) = [];
length_s(isCensored) = [];

durations = length_s(:);
starts = starts(:);
ends = ends(:);

% Assign Origin Phase based on the mask value at the bout's starting index
if isempty(starts)
    origins = strings(0, 1);
else
    start_mask_vals = mask(starts);
    origins = strings(size(starts));

    if mode == 2 % Light/Dark
        origins(start_mask_vals == 1) = "Light";
        origins(start_mask_vals == 0) = "Dark";
    elseif mode == 3 % Constant Darkness
        origins(start_mask_vals == 1) = "SubjDay";
        origins(start_mask_vals == 0) = "SubjNight";
    else % Unsegmented Entire Trace
        origins(:) = "Unsegmented";
    end
end
end

%% ------------------------------------------------------------------------
function mustBeActivityWave(x)
%MUSTBEACTIVITYWAVE Double vector of 0s and 1s, 0 = sleep, 1 = wake.
if ~isa(x, 'double')
    error('segment_sleep_phases:NotActivityWave', ...
        ['Input must be a DOUBLE activity wave (0 = sleep, 1 = wake), as ' ...
        'returned by sleeptrace; got %s. Logical input is refused because ' ...
        'the earlier convention here was logical true = SLEEP, and a stale ' ...
        'mask would be segmented with sleep and wake swapped, without ' ...
        'error. If your data are already 0 = sleep, 1 = wake, pass double(x).'], ...
        class(x));
end
if any(x ~= 0 & x ~= 1)          % NaN also fails this test
    error('segment_sleep_phases:NotBinary', ...
        'Activity wave must contain only 0 (sleep) and 1 (wake).');
end
end