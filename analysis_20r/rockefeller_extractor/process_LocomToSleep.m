function sleep_out = process_LocomToSleep(options)
% PROCESS_LOCOMTOSLEEP Converts raw locomotor time-series into segmented sleep trace data.
%
%   This function serves as the primary pipeline for processing Drosophila
%   locomotor activity data. It handles file selection, initial entrainment
%   trimming, automated death detection, sleep-wake conversion, and circadian
%   phase segmentation.
%
%   Pipeline Workflow:
%   1. File Selection: Prompts for a directory (skipped when a valid 'folderPath'
%      is supplied) and then for specific .txt files (skipped when 'selectAll'
%      is true). Set both to run with no dialogs at all.
%   2. Input Validation: Each file must read as a single row or column of finite
%      numbers. Empty files, multi-column exports (e.g. a whole DAM monitor
%      board, which would otherwise be flattened column-major into one chimeric
%      "fly") and files containing NaN/Inf are discarded with a warning.
%   3. Entrainment Trimming: Discards the first 'Teliminate' hours of data to
%      account for environmental acclimation.
%   4. Death Detection: Scans for 'deadThresh' consecutive hours of total
%      inactivity. If found, the data is truncated starting 24 hours *prior*
%      to the onset of death to exclude dying artifacts.
%   5. Length Validation: Discards flies that do not possess at least 'minDays'
%      of valid data after trimming and death truncation.
%   6. Sleep Conversion: Calls sleeptrace() to convert locomotor activity into
%      binary sleep/wake arrays (defined as 5+ minutes of inactivity).
%      sleeptrace returns an ACTIVITY wave, DOUBLE, 0 = sleep and 1 = wake.
%      That single convention is used unchanged by every step that follows --
%      ironout, segment_sleep_phases, and the returned masks -- so no inversion
%      happens anywhere in the pipeline. ironout and segment_sleep_phases refuse
%      logical input, because an earlier version of this code used logical
%      true = SLEEP and a stale mask of that kind would otherwise pass through
%      with sleep and wake silently swapped.
%   7. Brief-Wake Bridging: Calls ironout() to absorb wake runs no longer than
%      'bridgeWakeSeconds' that have sleep on both sides. ON BY DEFAULT; see the
%      WARNING below.
%   8. Phase Segmentation: Calls segment_sleep_phases() to divide the data into
%      circadian phases (e.g., Light/Dark) and calculate bout architecture.
%
%   ======================================================================
%   WARNING: BRIEF-WAKE BRIDGING IS ON BY DEFAULT AND CHANGES SLEEP SCORING
%   ======================================================================
%   By default ('bridgeWakeSeconds' = 1) any wake run lasting no more than one
%   second, with sleep on BOTH sides, is relabelled as sleep and the flanking
%   bouts are merged into one. The biological argument is that an isolated beam
%   break inside a long quiescent stretch is more likely a postural adjustment
%   than a genuine arousal.
%
%   This is a deliberate DEPARTURE from the conventional five-minute inactivity
%   definition implemented by sleeptrace. Samples in which the animal actually
%   moved are being scored as sleep. Consequences:
%
%     * Total sleep increases and wake decreases, slightly.
%     * Sleep bouts get longer, and there are fewer of them, because bouts that
%       were separate are now merged.
%     * THE WAKE-BOUT DISTRIBUTION IS TRUNCATED AT ITS SHORT END. Every
%       interior 1 s wake bout is removed by construction, so the shortest
%       wake bouts no longer exist in the output. Anything fitted to wake-bout
%       durations - exponential, power-law, hyperexponential, HMM dwell times -
%       is affected, and fits are NOT comparable to fits from unbridged data.
%       If you are modeling bout-length distributions, consider setting
%       bridgeWakeSeconds = 0 and handling arousals in the model instead.
%     * Results are not directly comparable to published work that uses the
%       standard five-minute rule without bridging.
%
%   MEASURED EFFECT on 1 s Drosophila DD recordings (median per fly):
%
%                              Canton-S (n=16)    per0 (n=148)
%       1 s wake bouts             17.9%             14.3%     of all wake bouts
%       wake bouts absorbed          50                44
%       total sleep change         +0.01 pp          +0.01 pp
%       wake bout count            -17.9%            -14.3%
%       sleep bout count           -17.9%            -14.3%
%       median sleep bout          +19.9%            +18.1%
%
%   Note how lopsided this is. Total sleep barely moves - a hundredth of a
%   percentage point - because only isolated seconds are relabelled. But nearly
%   one wake bout in six disappears and the median sleep bout grows by about a
%   fifth. Summary measures of sleep AMOUNT are essentially untouched; measures
%   of sleep ARCHITECTURE are not.
%
%   Set bridgeWakeSeconds = 0 to disable it and recover the conventional
%   definition exactly. The value used is recorded in
%   sleep_out.bridgeWakeSeconds, and the amount of relabelling per animal is
%   reported as flyN.bridgedWakeBouts and flyN.bridgedWakeSamples so the size
%   of the effect can be measured rather than assumed.
%
%   Bridging is applied to the activity wave only, AFTER sleeptrace and BEFORE
%   segmentation. It never alters the raw activity, and it cannot create a new
%   sleep bout: both flanks are already sleep, so an existing bout can only be
%   extended.
%
%   Error Isolation:
%   Every per-file problem raises a warning and skips to the next file, so one
%   unreadable or malformed record can never abort a batch. A fly is dropped,
%   never silently altered. Each warning names the fly by its position in the
%   SELECTED file list, e.g. 'Fly 19 ("monitor3_ch19.txt") discarded: ...', and
%   a roll-up of every discarded fly is printed at the end of the run. Note that
%   'Fly 19' is the 19th file you selected, which is NOT in general the same as
%   the flyN field in the output: rejected files consume no flyN index.
%
%   Discard warning identifiers, all suppressible with warning('off', id):
%       PIPELINE:EmptyFile            no numeric data could be read
%       PIPELINE:MultiColumnFile      not a single row/column of counts
%       PIPELINE:NonFiniteData        contains NaN or Inf
%       PIPELINE:ShorterThanTeliminate  file does not outlast the trim window
%       PIPELINE:DiedTooEarly         no complete 24 h pre-death window
%       PIPELINE:InsufficientData     less than minDays after all trimming
%       PIPELINE:SleepTraceFailed     sleeptrace declined to classify
%       PIPELINE:MissingDayLength     no entry in a per-fly dayLength map
%       PIPELINE:BadDayLength         per-fly entry is not a positive finite scalar
%       PIPELINE:FileFailed           any other unexpected error
%
%   Batch-level notices, raised once per call rather than per fly:
%       PIPELINE:PhaseAlignment       retained data does not begin at ZT0
%       PIPELINE:PerFlyDayLength      per-fly mode is active; daytime:nighttime
%                                     is being read as a ratio
%       PIPELINE:BridgeBelowResolution
%                                     bridgeWakeSeconds is shorter than one
%                                     sample, so no bridging can occur
%
%   One warning does NOT discard the fly:
%       PIPELINE:EmptyPhase           a phase received no bins at all, so its
%                                     mean() is NaN. Two routes reach this:
%
%                                     (a) A degenerate photoperiod: daytime=0 or
%                                     daytime=dayLength empties one phase at ANY
%                                     trace length. This is the route you are
%                                     likely to meet in practice, e.g. setting
%                                     daytime=0 for a DD run while leaving
%                                     segmentMode=2.
%
%                                     (b) A retained trace short enough to sit
%                                     entirely inside one phase. This requires
%                                     minDays*dayLength <= the longer phase
%                                     (minDays <= 0.5 under 12:12), because any
%                                     trace spanning a whole cycle necessarily
%                                     contains both phases. At the default
%                                     minDays=3 this route is unreachable: the
%                                     fly is dropped by PIPELINE:InsufficientData
%                                     before segmentation ever runs.
%
%   Syntax:
%       sleep_out = process_LocomToSleep() ; % 24 h, dialogs
%       sleep_out = process_LocomToSleep(Name=Value, ...)
%
%
%       D = 'C:\Users\sheyu\Box\Data from Shreenath\Selected_flies\CS_DD';
%       S = process_LocomToSleep(timeInterval=1, segmentMode=3, ...
%           folderPath=D, selectAll=true);    24 h, no dialogs
%
%       S = process_LocomToSleep(dayLength=18, daytime=9, nighttime=9, ...);  % any scalar
%
% *** If first screening for periodicity and then passsing that on to
% *** decide which flies' locomotion to convert to sleep
%       % D = 'C:\Users\sheyu\Box\Data from Shreenath\Selected_flies\CS_DD';
%       % f = dir(fullfile(D,'*.txt'));
%       % names = strings(numel(f),1);  tau = nan(numel(f),1);  pv = nan(numel(f),1);
%
% for k = 1:numel(f)
%     names(k) = string(f(k).name);
%     a = readmatrix(fullfile(D,f(k).name));  a = a(:);
%     a = a(24*3600+1:end);                       % match Teliminate
%     r = chi2p(a, InputInterval=seconds(1), BinInterval=minutes(6), ...
%         PeriodRange=hours([16 32]), BlockLength=hours(2), ...
%         NumSurrogates=3000);
%     tau(k) = r.TauHours;  pv(k) = r.PValue;
% end
%
% keep = pv <= 0.05;                              % your own rhythmicity criterion
% S = process_LocomToSleep(timeInterval=1, segmentMode=3, Teliminate=24, ...
%     dayLength = dictionary(names(keep), tau(keep)), ...
%     folderPath=D, selectAll=true);

%% *** using the output
%   sleepFraction = 1 - mean(f.entire);      % or nnz(f.entire==0)/numel(f.entire)
%   sleepMinPerDay = sleepFraction * 1440;
%   sleepBoutLengths = CalcBoutSizes(f.entire, 0);   % runs of ZEROS
%%


%   Name-Value Arguments:
%       timeInterval - (numeric) Duration of each sample bin in seconds. (Default: 1)
%       Teliminate   - (numeric) Hours of data to discard at the start.
%                      A non-multiple of 'dayLength' is permitted, not an error:
%                      the phase mask is realigned to
%                      ZT mod(startZT + Teliminate, dayLength) and a
%                      PIPELINE:PhaseAlignment warning records the new start
%                      phase. Choose it so that startZT + Teliminate is a
%                      multiple of dayLength if you want the retained data to
%                      begin at ZT0. (Default: 24)
%       deadThresh   - (numeric) Consecutive hours of zero activity required to
%                      flag a fly as dead. (Default: 12)
%                      Note: the amount of data discarded before the onset of
%                      death is fixed at 24 real hours, independent of dayLength.
%       minDays      - (numeric) Minimum days of valid data required post-trimming
%                      to include the fly in the final output. (Default: 3)
%
%                      A "day" here means one 'dayLength' cycle, NOT 24 absolute
%                      hours. This is intentional: the threshold scales with the
%                      experiment's own cycle, so minDays=3 requires 3 complete
%                      cycles whatever dayLength is. Under an 18 hr cycle that is
%                      54 hrs, not 72. Fractional values are allowed, and are the
%                      only way to retain windows shorter than one full cycle.
%       dayLength    - Total hours in one complete circadian cycle. Either a
%                      positive scalar applied to every fly (Default: 24), or a
%                      PER-FLY map giving each animal its own period.
%
%                      Per-fly forms, all keyed by file name:
%                          dictionary(["a.txt" "b.txt"], [23.8 24.3])
%                          containers.Map({'a.txt','b.txt'}, {23.8, 24.3})
%                          table(fileNames, tauHours)
%
%                      Intended for DD, where each fly free-runs at its own tau.
%                      A single 24 h mask drifts out of register with a fly whose
%                      tau is not 24: at tau = 24.6 over 9 days the mask ends up
%                      about 4.5 h out of phase.
%
%                      HOW MUCH THIS ACTUALLY BUYS YOU. Measured on 14 rhythmic
%                      Canton-S DD records of about 9 cycles, with the phase
%                      offset optimised separately so period and phase are not
%                      confounded, the subjective day-night sleep amplitude was
%                      0.282 at a fixed 24 h, 0.279 at each fly's measured tau,
%                      and 0.287 at the best-fitting tau. That is a 2% gain,
%                      while allowing tau to float freely inflates the same
%                      statistic by 58% on ARRHYTHMIC per0 records, which is
%                      pure overfitting. Over records this short the amplitude
%                      surface is nearly flat in tau, so per-fly periods are not
%                      where the error lives.
%
%                      The dominant error is startZT, not the period. On those
%                      same records the amplitude roughly doubled (0.14 to 0.28)
%                      when the phase offset was fitted rather than assumed,
%                      meaning the default startZT=0 was substantially wrong for
%                      them. Establish startZT before reaching for per-fly tau.
%
%                      Per-fly tau earns its keep when the period departs
%                      further from nominal or the record runs longer, so that
%                      accumulated drift becomes large relative to a phase
%                      window: short-period and long-period mutants, or
%                      multi-week recordings.
%
%                      Three things change in per-fly mode:
%                        * daytime and nighttime are read as a RATIO, not as
%                          absolute hours, and rescaled to each fly's cycle.
%                          12:12 means "the first half of this fly's own cycle",
%                          which is the circadian-time convention. Their sum need
%                          not equal anything.
%                        * minDays counts cycles of that fly's tau, so the
%                          length threshold differs slightly between flies.
%                        * the phase offset is mod(startZT + Teliminate, tau).
%
%                      A file absent from the map is DISCARDED with
%                      PIPELINE:MissingDayLength rather than analysed at some
%                      default. If a fly's period could not be measured, no phase
%                      assignment for it is trustworthy; include an explicit
%                      entry (e.g. 24) to analyse it at a nominal cycle anyway.
%
%                      Batch-level metadata that has no single correct value in
%                      this mode (.onewholeday, .dataStartZT) is set to NaN. The
%                      values actually used are recorded per fly as
%                      .dayLengthUsed, .daytimeUsed and .dataStartZT.
%       daytime      - (numeric) Duration of the light (or subjective day) phase
%                      in hours. (Default: 12)
%       nighttime    - (numeric) Duration of the dark (or subjective night) phase
%                      in hours. Must satisfy daytime + nighttime == dayLength;
%                      an inconsistent triple is now rejected up front. (Default: 12)
%       segmentMode  - (numeric) How to segment the phases downstream.
%                      1: Entire trace only.
%                      2: Light / Dark split. (Default)
%                      3: Subjective Day / Subjective Night (DD).
%       startZT      - (numeric) The ZT/CT hour at which the RAW recording's first
%                      sample was taken, before any trimming. (Default: 0)
%
%                      This is the pipeline's only assumption about absolute
%                      time, and it cannot be checked: activity counts carry no
%                      clock. Together with 'Teliminate' it fixes the starting
%                      phase of the retained data,
%                      ZT mod(startZT + Teliminate, dayLength), which is passed
%                      to segment_sleep_phases as 'phaseOffset' and recorded in
%                      the output as '.dataStartZT'.
%
%                      Getting it wrong does not raise an error. It shifts every
%                      Light/Dark assignment, and every bout's 'origin_phase', by
%                      a constant, while leaving the totals and the trace itself
%                      untouched. The default of 0 assumes the recording began at
%                      lights-on; set it explicitly if it did not.
%       folderPath   - (char/string) Path to the data directory. If empty, a UI
%                      folder selection dialog will open. (Default: '')
%       bridgeWakeSeconds
%                    - (numeric) Maximum duration of a wake run that will be
%                      absorbed into sleep when it has sleep on both sides.
%                      (Default: 1 second. Set 0 to disable.)
%
%                      See the WARNING above before changing or relying on this.
%                      Only INTERIOR runs are bridged: a short wake run touching
%                      the first or last sample has sleep on one side only, and
%                      its true extent is censored by the recording boundary.
%
%                      The threshold is in SECONDS and is converted with
%                      floor(bridgeWakeSeconds / timeInterval), so coarse
%                      sampling disables it automatically - at timeInterval = 60
%                      no run can be 1 s long. A PIPELINE:BridgeBelowResolution
%                      warning is raised in that case rather than failing quietly.
%       selectAll    - (logical) Skip the file-selection dialog and process every
%                      .txt file in the folder. Set true together with folderPath
%                      to run the whole pipeline non-interactively, e.g. from a
%                      batch script or an automated test. (Default: false)
%
%   Output Structure (sleep_out):
%       .onewholeday   - Metadata: 'dayLength' value used.
%       .daytime       - Metadata: 'daytime' value used.
%       .nighttime     - Metadata: 'nighttime' value used.
%       .Teliminate    - Metadata: 'Teliminate' value used.
%       .timeInterval  - Metadata: 'timeInterval' value used.
%       .segmentMode   - Metadata: 'segmentMode' value used.
%       .startZT       - Metadata: 'startZT' value used (raw recording start).
%       .dataStartZT   - Metadata: the ZT/CT hour at which every fly's '.entire'
%                        trace begins, i.e. mod(startZT + Teliminate, dayLength).
%                        Records the phase alignment a saved structure was built
%                        under, so results remain interpretable later.
%                        NaN in per-fly dayLength mode; see flyN.dataStartZT.
%       .outputCoding  - Metadata: always "sleepIs0" (0 = sleep, 1 = wake,
%                        DOUBLE). Kept so saved structures, and any code that
%                        reads this field, remain self-describing.
%       .bridgeWakeSeconds
%                      - Metadata: brief-wake bridging threshold in seconds that
%                        was applied (0 means none). See the WARNING above.
%       .perFlyDayLength
%                      - Metadata: true when a per-fly dayLength map was used.
%       .daytimeFraction
%                      - Metadata: fraction of each cycle treated as subjective
%                        day. Present in per-fly mode only, where daytime and
%                        nighttime are read as a ratio.
%       .includedFiles - Cell array of filenames that successfully passed all filters.
%       .fly1, .fly2...- Nested structures for each valid fly containing the binary
%                        masks (e.g., .entire, .light, .dark) and bout architecture
%                        stats (.stats.sleep_bouts, .stats.wake_bouts).
%                        The masks are DOUBLE with 0 = sleep and 1 = wake, the
%                        sleeptrace convention, so mean(.entire) is the WAKE
%                        fraction and 1 - mean(.entire) the sleep fraction.
%                        This matches the older scripts in this repo
%                        (boutsONOFF, BoutDistribs), so the masks drop straight
%                        into them. Note that CalcBoutSizes REQUIRES double: it
%                        assigns sentinel values (100/101) into the array, which
%                        collapse to true/false on a logical input.
%                        .stats is named semantically: .stats.sleep_bouts are
%                        runs of 0s, .stats.wake_bouts runs of 1s, and their
%                        .starts/.ends index directly into .entire.
%
%                        Sanity check: every interior sleep bout is at least
%                        5 minutes long by definition, while wake bouts have no
%                        minimum. Short "sleep" runs mean the coding is not
%                        what you assumed.
%                        flyN and includedFiles{N} always refer to the same file;
%                        rejected files consume neither index.
%
%                        Each flyN also carries the settings it was analysed
%                        under, so a saved structure stands alone:
%                          .dayLengthUsed       cycle length for this animal
%                          .daytimeUsed         subjective day length, hours
%                          .dataStartZT         ZT/CT at which .entire begins
%                          .bridgedWakeBouts    wake runs absorbed by ironout
%                          .bridgedWakeSamples  samples that amounted to
%
%   Dependencies:
%       sleeptrace.m
%       ironout.m               (brief-wake bridging; see the WARNING above)
%       segment_sleep_phases.m
%
%   Example 1: Standard 24h LD cycle (default parameters)
%       sleep_data = process_LocomToSleep();
%
%   Example 2: 9:9 LD cycle with 60-second bins, skipping the UI folder dialog
%       sleep_data = process_LocomToSleep(timeInterval=60, dayLength=18, ...
%                                         daytime=9, nighttime=9, ...
%                                         folderPath='C:\Data\FlyRun1');
%
%   Example 3: Fully non-interactive, for a batch script or an automated test
%       sleep_data = process_LocomToSleep(timeInterval=60, ...
%                                         folderPath='C:\Data\FlyRun1', ...
%                                         selectAll=true);
%
%   Example 4: Total sleep per phase, in minutes, for every surviving fly
%       for k = 1:numel(sleep_data.includedFiles)
%           f = sleep_data.(sprintf('fly%d', k));
%           mins = sleep_data.timeInterval / 60;
%           fprintf('%-24s light %7.1f min   dark %7.1f min\n', ...
%               sleep_data.includedFiles{k}, sum(f.light)*mins, sum(f.dark)*mins);
%       end
%
%   Example 5: Sleep FRACTION per phase. Use mean(), not sum() with a
%   hand-built denominator: numel(f.light) and numel(f.dark) are generally
%   NOT equal, because the trace rarely spans a whole number of cycles and
%   because a non-multiple 'Teliminate' shifts the starting phase.
%       f = sleep_data.fly1;
%       frac_light = mean(f.light);     % divides by numel(f.light) for you
%       frac_dark  = mean(f.dark);
%       fprintf('light %.3f (%d bins)   dark %.3f (%d bins)\n', ...
%           frac_light, numel(f.light), frac_dark, numel(f.dark));
%
%   Errors:
%       PIPELINE:PhaseDurationMismatch  daytime + nighttime ~= dayLength
%                                       (scalar dayLength mode only)
%       PIPELINE:PhaseRatioUndefined    daytime + nighttime is not positive, so
%                                       the day:night ratio is undefined
%                                       (per-fly dayLength mode only)
%       PIPELINE:BadDayLength           dayLength is neither a positive scalar
%                                       nor a valid per-fly map
%       Both remaining hard errors concern folder/file selection: an empty
%       directory, or a cancelled dialog. Everything else warns and skips.

arguments
    options.timeInterval (1,1) {mustBeNumeric, mustBePositive} = 1
    options.Teliminate   (1,1) {mustBeNumeric, mustBeNonnegative} = 24
    options.deadThresh   (1,1) {mustBeNumeric, mustBePositive} = 12
    options.minDays      (1,1) {mustBeNumeric, mustBePositive} = 3
    options.dayLength = 24   % scalar, or a per-fly map: see localResolveDayLength
    options.daytime      (1,1) {mustBeNumeric, mustBeNonnegative} = 12
    options.nighttime    (1,1) {mustBeNumeric, mustBeNonnegative} = 12
    options.segmentMode  (1,1) {mustBeMember(options.segmentMode, [1, 2, 3])} = 2
    options.startZT      (1,1) {mustBeNumeric, mustBeNonnegative} = 0
    options.folderPath   (1,:) char = ''
    options.selectAll    (1,1) logical = false
    options.bridgeWakeSeconds (1,1) {mustBeNumeric, mustBeNonnegative} = 1
end

% --- Resolve dayLength: one global cycle, or one per fly ---
[perFlyTau, tauLookup, scalarDayLength] = localResolveDayLength(options.dayLength);

if perFlyTau
    % Per-fly mode. Each fly gets its own free-running period, so the
    % light/dark split cannot be given in absolute hours: 12 h of a 23.6 h
    % cycle and 12 h of a 24.6 h cycle are different fractions of their
    % cycles. The RATIO daytime:nighttime is used instead and rescaled to
    % each fly's tau, which is the standard circadian-time convention
    % (subjective day is the first half of the animal's own cycle).
    if (options.daytime + options.nighttime) <= 0
        error('PIPELINE:PhaseRatioUndefined', ...
            ['In per-fly dayLength mode the daytime:nighttime RATIO defines ' ...
            'the split, so their sum must be positive.']);
    end
    daytimeFraction = options.daytime / (options.daytime + options.nighttime);
    warning('PIPELINE:PerFlyDayLength', ...
        ['Per-fly dayLength supplied for %d file(s). daytime:nighttime is ' ...
        'being read as a ratio (%.3f of each cycle is subjective day) and ' ...
        'rescaled to each fly''s own period. minDays and the phase offset ' ...
        'scale with it too; see .dayLengthUsed on each fly.'], ...
        numEntries(tauLookup), daytimeFraction);
else
    % Scalar mode, unchanged: daytime and nighttime must tile the cycle.
    daytimeFraction = NaN;
    if abs((options.daytime + options.nighttime) - scalarDayLength) > 1e-9
        error('PIPELINE:PhaseDurationMismatch', ...
            ['daytime (%g hrs) + nighttime (%g hrs) = %g hrs, which does not match ' ...
            'dayLength (%g hrs). These must tile exactly one circadian cycle.'], ...
            options.daytime, options.nighttime, ...
            options.daytime + options.nighttime, scalarDayLength);
    end

    % --- Starting Phase of the Retained Data ---
    % The raw recording begins at ZT/CT 'startZT'; 'Teliminate' hours are then
    % discarded, so the retained data begins this far into the cycle. This single
    % value carries the pipeline's only untestable assumption about absolute time:
    % nothing in the activity counts reveals when the recording started, so a
    % wrong startZT shifts every phase assignment uniformly and silently.
    dataStartZT = mod(options.startZT + options.Teliminate, scalarDayLength);

    % --- Phase Alignment Notice ---
    % Not an error: the mask is realigned via phaseOffset below. This only flags
    % that the retained data no longer begins at ZT0.
    if dataStartZT ~= 0
        warning('PIPELINE:PhaseAlignment', ...
            ['Recording starts at ZT%g and Teliminate discards %g hrs, so the ' ...
            'retained data begins at ZT%g of a %g hr cycle. The phase mask has ' ...
            'been shifted accordingly.'], ...
            options.startZT, options.Teliminate, dataStartZT, scalarDayLength);
    end
end

% 1. Folder Selection UI
folderPath = options.folderPath;
if isempty(folderPath) || ~isfolder(folderPath)
    folderPath = uigetdir('', 'Select Directory Containing Fly Data Files');
    if isequal(folderPath, 0)
        error('Process canceled: No directory selected.');
    end
end

txtFiles = dir(fullfile(folderPath, '*.txt'));
if isempty(txtFiles)
    error('No .txt files found in folder: %s', folderPath);
end

% 2. Select Files (via dialog, unless selectAll was requested)
if options.selectAll
    selectedFiles = txtFiles;
else
    fileNames = {txtFiles.name};
    [selectedIdx, tf] = listdlg(...
        'PromptString', 'Select fly data files to include:', ...
        'SelectionMode', 'multiple', ...
        'ListSize', [300, 400], ...
        'ListString', fileNames);

    if ~tf || isempty(selectedIdx)
        error('Process canceled: No files selected.');
    end
    selectedFiles = txtFiles(selectedIdx);
end

% Conversion constants (hours -> sample bins)
bins_per_hour    = 3600 / options.timeInterval;
Telim_bins       = round(options.Teliminate * bins_per_hour);
dead_thresh_bins = round(options.deadThresh * bins_per_hour);
pre_death_bins   = round(24 * bins_per_hour);
% A wake run qualifies for bridging when its duration does not exceed
% bridgeWakeSeconds, i.e. when it spans at most floor(seconds/timeInterval)
% samples. Coarse sampling therefore disables it automatically: at
% timeInterval = 60 s no run can be 1 s long, so bridge_bins is 0.
bridge_bins = floor(options.bridgeWakeSeconds / options.timeInterval);
if options.bridgeWakeSeconds > 0 && bridge_bins < 1
    warning('PIPELINE:BridgeBelowResolution', ...
        ['bridgeWakeSeconds (%g s) is shorter than one sample (%g s), so no ' ...
        'wake bout can qualify and no bridging will occur.'], ...
        options.bridgeWakeSeconds, options.timeInterval);
end

% min_req_bins depends on the cycle length, so in per-fly mode it is
% resolved inside the loop instead.
if ~perFlyTau
    min_req_bins = round(options.minDays * scalarDayLength * bins_per_hour);
end

% Store output structure metadata. In per-fly mode the batch-level cycle
% fields are set to NaN on purpose: there is no single correct value, and
% leaving a plausible-looking scalar there would invite downstream code to
% use the wrong cycle silently. The real values live on each fly as
% .dayLengthUsed / .dataStartZT.
sleep_out = struct();
if perFlyTau
    sleep_out.onewholeday   = NaN;
    sleep_out.daytimeFraction = daytimeFraction;
    sleep_out.perFlyDayLength = true;
else
    sleep_out.onewholeday   = scalarDayLength;
    sleep_out.perFlyDayLength = false;
end
sleep_out.daytime      = options.daytime;
sleep_out.nighttime    = options.nighttime;
sleep_out.Teliminate   = options.Teliminate;
sleep_out.timeInterval = options.timeInterval;
sleep_out.segmentMode  = options.segmentMode;
sleep_out.startZT      = options.startZT;
% Recorded so a saved structure is self-describing: nothing else in the
% output reveals which way round the binary masks are. Fixed: the whole
% pipeline uses the sleeptrace convention, 0 = sleep, 1 = wake.
sleep_out.outputCoding = "sleepIs0";
sleep_out.bridgeWakeSeconds = options.bridgeWakeSeconds;
if perFlyTau
    sleep_out.dataStartZT = NaN;        % per fly; see flyN.dataStartZT
else
    sleep_out.dataStartZT = dataStartZT; % ZT at which every .entire begins
end
sleep_out.includedFiles = {};

validFlyCount = 0;
keptIdx = [];   % indices into selectedFiles that survived every filter

% 3. File Processing Loop
% Every rejection path warns and skips to the next file, so one unreadable or
% malformed fly can never abort the whole batch. Warnings identify the fly by
% its position in the SELECTED file list, so 'Fly 19' is the 19th entry you
% picked in the dialog (or the 19th .txt in the folder when selectAll=true),
% which is not in general the same as its flyN field in the output.
for i = 1:length(selectedFiles)
    fileName = selectedFiles(i).name;
    flyTag   = sprintf('Fly %d ("%s")', i, fileName);
    try
        filePath = fullfile(selectedFiles(i).folder, fileName);
        raw_activity = readmatrix(filePath);

        % Step 0: Validate the shape and content of the file
        if isempty(raw_activity)
            warning('PIPELINE:EmptyFile', ...
                '%s discarded: no numeric data could be read.', flyTag);
            continue;
        end
        if ~isvector(raw_activity)
            % A multi-column export (e.g. a whole DAM monitor board) would be
            % flattened column-major by (:), silently splicing several channels
            % into one "fly". Reject it rather than analyse a chimera.
            warning('PIPELINE:MultiColumnFile', ...
                ['%s discarded: read as a %s array. Expected a single row ' ...
                'or column of activity counts for one fly.'], ...
                flyTag, mat2str(size(raw_activity)));
            continue;
        end
        raw_activity = raw_activity(:); % Force column vector

        % Resolve this fly's cycle. A fly absent from a per-fly map is DROPPED
        % rather than quietly analysed at some default: if its period could not
        % be measured, no phase assignment for it is trustworthy, and silently
        % substituting 24 h is exactly the hidden assumption this mode exists to
        % remove. Include an explicit entry to analyse a fly at a nominal cycle.
        if perFlyTau
            if ~isKey(tauLookup, string(fileName))
                warning('PIPELINE:MissingDayLength', ...
                    ['%s discarded: no per-fly dayLength supplied. Add an entry ' ...
                    '(e.g. 24) to analyse it at a nominal cycle.'], flyTag);
                continue;
            end
            flyDayLength = double(tauLookup(string(fileName)));
            if ~isscalar(flyDayLength) || ~isfinite(flyDayLength) || flyDayLength <= 0
                warning('PIPELINE:BadDayLength', ...
                    '%s discarded: per-fly dayLength must be a positive finite scalar.', ...
                    flyTag);
                continue;
            end
            flyDaytime   = daytimeFraction * flyDayLength;
            flyStartPhase = mod(options.startZT + options.Teliminate, flyDayLength);
            min_req_bins = round(options.minDays * flyDayLength * bins_per_hour);
        else
            flyDayLength  = scalarDayLength;
            flyDaytime    = options.daytime;
            flyStartPhase = dataStartZT;
        end

        if any(~isfinite(raw_activity))
            warning('PIPELINE:NonFiniteData', ...
                ['%s discarded: contains %d non-finite value(s) (NaN/Inf), ' ...
                'usually unparsed text or blank fields.'], ...
                flyTag, nnz(~isfinite(raw_activity)));
            continue;
        end

        % Step A: Cut initial entrainment window
        if length(raw_activity) <= Telim_bins
            warning('PIPELINE:ShorterThanTeliminate', ...
                ['%s discarded: only %.1f hrs long, which does not exceed the ' ...
                '%.1f hr Teliminate window.'], ...
                flyTag, length(raw_activity)/bins_per_hour, options.Teliminate);
            continue;
        end
        activity = raw_activity(Telim_bins + 1 : end);

        % Step B: Death detection & truncation (24h prior to death start)
        inactive = (activity <= 0);
        transitions = diff([false; inactive; false]);
        boutStarts  = find(transitions == 1);
        boutEnds    = find(transitions == -1) - 1;
        boutDurations = boutEnds - boutStarts + 1;

        deathBoutIdx = find(boutDurations >= dead_thresh_bins, 1, 'first');

        if ~isempty(deathBoutIdx)
            death_start_sample = boutStarts(deathBoutIdx);
            cut_sample = death_start_sample - pre_death_bins;

            if cut_sample <= 1
                warning('PIPELINE:DiedTooEarly', ...
                    ['%s discarded: death onset at %.1f hrs post-trim leaves no ' ...
                    'complete 24 hr pre-death window to retain.'], ...
                    flyTag, death_start_sample/bins_per_hour);
                continue;
            else
                activity = activity(1 : cut_sample - 1);
            end
        end

        % Step C: Minimum length validation
        if length(activity) < min_req_bins
            warning('PIPELINE:InsufficientData', ...
                ['%s discarded: remaining data (%.1f hrs) is less than the ' ...
                'required %g days (%.1f hrs).'], ...
                flyTag, length(activity)/bins_per_hour, options.minDays, ...
                options.minDays * flyDayLength);
            continue;
        end

        % Step D: Sleep-wake conversion.
        % sleeptrace returns an activity wave, DOUBLE, 0 = sleep, 1 = wake. This
        % is the convention for the entire pipeline; it is passed on unchanged.
        activity_wave = sleeptrace(activity, options.timeInterval);
        if ~isequal(size(activity_wave), size(activity)) || any(isnan(activity_wave))
            % sleeptrace signals a refusal to classify by returning a scalar NaN.
            warning('PIPELINE:SleepTraceFailed', ...
                ['%s discarded: sleeptrace declined to classify this trace ' ...
                '(most often negative activity values).'], flyTag);
            continue;
        end

        % Step D2: Bridge brief wake interruptions. A wake run no longer than
        % 'bridgeWakeSeconds', flanked by sleep on BOTH sides, is relabelled as
        % sleep, merging the two flanking bouts into one. The rationale is that a
        % single beam break inside a long quiescent stretch is more likely a
        % postural shift than an arousal. See the WARNING in the header: this
        % departs from the conventional 5-minute definition and reshapes the
        % wake-bout distribution.
        % Applied here, BEFORE segmentation: the whole point is to merge bouts,
        % so the phase masks and .stats must be computed from the bridged trace.
        [activity_wave, nBridged, nBridgedSamples] = ...
            ironout(activity_wave, options.bridgeWakeSeconds, options.timeInterval);

        % Step E: Circadian phase segmentation.
        % nighttime is only forwarded in scalar mode, where the caller supplied
        % it and a cross-check is meaningful. In per-fly mode the dark phase is
        % the exact remainder of that fly's cycle, so passing a separately
        % rounded value could trip the segmenter's consistency check by one bin.
        if perFlyTau
            flyNighttime = [];
        else
            flyNighttime = options.nighttime;
        end
        segmented_data = segment_sleep_phases(activity_wave, ...
            mode = options.segmentMode, ...
            timeInterval = options.timeInterval, ...
            dayLength = flyDayLength, ...
            daytime = flyDaytime, ...
            nighttime = flyNighttime, ...
            phaseOffset = flyStartPhase);

        % Record the cycle this fly was actually analysed under, so a saved
        % structure remains interpretable without the original call.
        segmented_data.dayLengthUsed = flyDayLength;
        segmented_data.daytimeUsed   = flyDaytime;
        segmented_data.dataStartZT   = flyStartPhase;
        % How much the bridging step actually changed this fly, so the effect can
        % be quantified rather than assumed.
        segmented_data.bridgedWakeBouts   = nBridged;
        segmented_data.bridgedWakeSamples = nBridgedSamples;

        % Step E2: Flag phases that received no bins at all. This is not an error
        % and the fly is kept, but mean() of an empty phase is NaN, which would
        % otherwise propagate silently into any group statistic.
        %
        % Reached either by a degenerate photoperiod (daytime=0 or
        % daytime=dayLength, which empties a phase at any length) or by a trace
        % short enough to sit inside one phase. The latter needs
        % minDays*dayLength <= the longer phase, since a trace spanning a whole
        % cycle always contains both; at the default minDays=3 the length filter
        % above fires first, so only the degenerate-photoperiod route survives.
        %
        % The check lives here rather than in segment_sleep_phases so the message
        % can name the fly, and so the segmenter stays silent by default.
        for fld = ["light" "dark" "subj_day" "subj_night"]
            if isfield(segmented_data, fld) && isempty(segmented_data.(fld))
                warning('PIPELINE:EmptyPhase', ...
                    ['%s: the ''%s'' phase contains no bins, so mean(.%s) is NaN. ' ...
                    'The retained %.1f hrs begin at ZT%g of a %g hr cycle and ' ...
                    'never enter that phase.'], ...
                    flyTag, fld, fld, numel(activity_wave)/bins_per_hour, ...
                    flyStartPhase, flyDayLength);
            end
        end

        % Step F: Store structure. Compute the index first so a failure here can
        % never leave a gap between flyN and includedFiles{N}.
        nextIdx = validFlyCount + 1;
        sleep_out.(sprintf('fly%d', nextIdx)) = segmented_data;
        sleep_out.includedFiles{nextIdx} = fileName;
        validFlyCount = nextIdx;
        keptIdx(end+1) = i; %#ok<AGROW>

    catch ME
        warning('PIPELINE:FileFailed', ...
            '%s discarded: unexpected error (%s): %s', ...
            flyTag, ME.identifier, ME.message);
    end
end

% 4. Roll-up. The per-fly warnings above say WHY each was dropped; this is a
% single place to see WHICH, without scrolling back through them.
nSelected = length(selectedFiles);
fprintf('Successfully processed %d out of %d selected files.\n', validFlyCount, nSelected);
discardedIdx = setdiff(1:nSelected, keptIdx);
if ~isempty(discardedIdx)
    fprintf('Discarded %d of %d:\n', numel(discardedIdx), nSelected);
    for k = discardedIdx
        fprintf('    Fly %d ("%s")\n', k, selectedFiles(k).name);
    end
end
end

%% ------------------------------------------------------------------------
function [perFly, lookup, scalarValue] = localResolveDayLength(dl)
%LOCALRESOLVEDAYLENGTH Accept a scalar cycle or a per-file map of cycles.
%
%   Returns perFly=false and a validated scalar, or perFly=true and a
%   dictionary keyed by file name (string) holding each fly's cycle in hours.
%   Accepted per-fly forms: dictionary, containers.Map, or a two-column table
%   whose first column is the file name and second the period in hours.

if isnumeric(dl)
    if ~isscalar(dl) || ~isfinite(dl) || dl <= 0
        error('PIPELINE:BadDayLength', ...
            ['Numeric dayLength must be a positive finite scalar. For per-fly ' ...
            'periods supply a dictionary, containers.Map, or two-column table ' ...
            'keyed by file name.']);
    end
    perFly = false; lookup = dictionary; scalarValue = double(dl);
    return
end

perFly = true; scalarValue = NaN;
if isa(dl, 'dictionary')
    lookup = dl;
elseif isa(dl, 'containers.Map')
    lookup = dictionary(string(keys(dl)), cell2mat(values(dl)));
elseif istable(dl)
    if width(dl) < 2
        error('PIPELINE:BadDayLength', ...
            'Table dayLength needs two columns: file name, then period in hours.');
    end
    lookup = dictionary(string(dl{:,1}), double(dl{:,2}));
else
    error('PIPELINE:BadDayLength', ...
        ['dayLength must be a positive scalar, or a dictionary, containers.Map ' ...
        'or two-column table keyed by file name. Got a %s.'], class(dl));
end

if numEntries(lookup) == 0
    error('PIPELINE:BadDayLength', 'Per-fly dayLength map is empty.');
end
v = values(lookup);
if ~isnumeric(v) || any(~isfinite(v)) || any(v <= 0)
    error('PIPELINE:BadDayLength', ...
        'Every per-fly dayLength must be a positive finite number of hours.');
end
end