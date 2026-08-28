function result = detect_shuffling(source, opts)
% DETECT_SHUFFLING Find and characterise shuffling episodes in one recording.
%
%   result = detect_shuffling(filepath)
%   result = detect_shuffling(filepath, opts)
%   result = detect_shuffling(preProcessV2Output, opts)
%   result = detect_shuffling(rawMatrix, opts)
%
%   runs the whole per-recording pass: clean the raw file, classify every
%   channel, detect shuffling episodes in the channels worth analysing,
%   promote extreme episodes to a hardware-artifact status, and optionally
%   score each episode against sleep/wake state.
%
%   Shuffling is sustained alternation between two adjacent beams
%   (15,14,15,14,...). It is detected on the cleaned position trace at
%   native resolution, upstream of multibeam_locomotion, because
%   sum(abs(diff(position))) cannot separate oscillating in place from
%   walking. The same detector does double duty as QC: near-continuous
%   alternation across a whole recording is a hardware fault, and the
%   existing death detection structurally cannot catch it, since that
%   looks for inactivity and this is continuous activity.
%
%   INPUTS
%   ------
%   source
%       Path to a raw .txt file, a PreProcessV2 output struct, or a raw
%       numeric matrix.
%
%   opts
%       Optional struct. All fields optional except as noted.
%
%       Acquisition
%         .N                   Channels to analyse. Default 64.
%         .Dt                  Seconds per sample. Default 0.2 (5 Hz).
%
%       Detection
%         .MinAlternations     Minimum runs per reported episode.
%                              Default 3. Detection is independent of
%                              this value except for the final gate, so
%                              set it once at the lowest threshold of
%                              interest and let shuffling_stats filter
%                              upward on NRuns; that is exactly
%                              equivalent to re-detecting.
%
%       QC (see shuffling_qc)
%         .Window              NChannels x 2 analysis window, raw file
%                              rows. STRONGLY RECOMMENDED -- see below.
%         .DiedMidRecording    NChannels x 1 logical from upstream.
%         .FlyLabels           1 x NChannels cell of names.
%
%       ChatteringArtifact thresholds
%         .ChatterMinCoreHours Single-episode core duration at or above
%                              which a fly is called a hardware artifact.
%                              Default 1 hour.
%         .ChatterMinFraction  Single-episode core duration as a fraction
%                              of the analysis window, at or above which
%                              the same applies. Default 0.20.
%
%       Sleep overlap (optional; omit to skip)
%         .Sleep               Struct, see "SLEEP OVERLAP" below. Every
%                              field of it is required when it is given.
%
%       Reporting
%         .OutputDir           Directory to write the QC and episode
%                              tables into as tab-delimited .txt. Default
%                              '' (do not write).
%         .Verbose             Print a summary. Default false.
%
%   OUTPUT
%   ------
%   result
%       Struct with fields:
%         .QC        Per-fly status table from shuffling_qc, with any
%                    ChatteringArtifact rows applied.
%         .Episodes  Struct of column vectors, one row per episode, with
%                    the fields listed under "EPISODE FIELDS".
%         .PerFly    Struct of column vectors, one row per fly:
%                    .NEpisodes, .ShuffleSamples (deduplicated),
%                    .ShuffleSeconds, .ShuffleFraction, .MaxNRuns,
%                    .MaxCoreSec, .WindowSeconds.
%         .Params    Every resolved parameter, for the record.
%         .Notes     Cell of strings describing what was and was not done.
%         .Sleep     The normalised sleep trace (polarity already
%                    resolved to a logical is-sleep vector per fly), or
%                    [] when none was supplied. Kept so that
%                    validate_shuffling can re-use it without the
%                    polarity having to be declared a second time.
%
%   EPISODE FIELDS
%   --------------
%     .FlyID                  Channel index
%     .FlyLabel               Channel name
%     .StartSample            First sample, as a RAW FILE ROW
%     .EndSample              Last sample, as a RAW FILE ROW
%     .StartSampleInWindow    First sample, relative to the fly's window
%     .StartTimeHours         Elapsed hours since the start of this fly's
%                             trimmed recording. Elapsed time, NOT ZT --
%                             circadian phase is deferred, see below.
%     .DurationSec            Outer duration, StartSample..EndSample
%     .CoreDurationSec        Alternation-only duration, excluding the
%                             dwell in the first and last runs. Use this
%                             one for statistics; see find_shuffles.
%     .NRuns                  Number of runs in the alternation
%     .BeamA, .BeamB          The two beams, in the order first visited.
%                             Not biologically meaningful, retained
%                             because a beam pair recurring across
%                             experiments is a hardware lead.
%     .MaxCoreRunSec          Longest interior run
%     .OverlapsSleep          True/false, or NaN when no sleep trace
%     .OverlapsWake           True/false, or NaN when no sleep trace.
%                             Both can be true: an episode may span a
%                             transition, and is not forced to one label.
%
%   THE WINDOW, AND WHY IT MATTERS MORE THAN ANYTHING ELSE HERE
%   -----------------------------------------------------------
%   process_LocomToSleep trims for entrainment and truncates on death
%   before it computes anything, so sleep_out sample indices do not
%   correspond 1:1 to raw file rows. Every index this function reports is
%   an explicit raw file row, and .StartTimeHours is measured from the
%   fly's own window start, so the two conventions cannot be confused by
%   accident. But that only holds if .Window is the same window
%   process_LocomToSleep used. Re-derive it from the options that call
%   was made with, or better, have that function hand it over directly.
%   Omitting .Window analyses whole files and flags every fly
%   WindowNotSupplied; a misalignment here would corrupt every overlap
%   statistic without raising an error.
%
%   SLEEP OVERLAP
%   -------------
%   opts.Sleep, when supplied, must have all four fields:
%
%     .Trace       nBins x NChannels numeric/logical matrix, or a
%                  1 x NChannels cell array of per-fly vectors (use the
%                  cell form when trimming left flies with different
%                  lengths).
%     .BinSeconds  Seconds per bin of that trace. Scalar, or one per fly.
%     .StartSample Raw file row that bin 1 corresponds to. Scalar, or one
%                  per fly. This is the trim offset, made explicit.
%     .Polarity    'ZeroIsSleep' or 'OneIsSleep'. REQUIRED, no default.
%
%   Polarity has no default on purpose. sleeptrace returns an activity
%   wave in which 0 means sleep, which is the opposite of what the name
%   suggests, and both its docstring and process_LocomToSleep's warn
%   about the confusion repeatedly. A wrong guess here would invert every
%   sleep/wake statistic silently, so the caller has to state it.
%
%   NOT DONE HERE
%   -------------
%   Day/night and any other ZT or circadian-phase split is deferred: this
%   pass is unphased, and reports elapsed hours since the start of each
%   fly's trimmed recording. Adding phase back needs a phase mask that
%   segment_sleep_phases builds internally but does not return.
%
%   See also FIND_SHUFFLES, SHUFFLING_QC, PREPROCESSV2, SHUFFLING_STATS.

    %% Resolve options

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    if ~isstruct(opts) || ~isscalar(opts)
        error('detect_shuffling:InvalidOpts', 'opts must be a scalar struct.');
    end

    N                  = getOpt(opts, 'N', 64);
    dt                 = getOpt(opts, 'Dt', 0.2);
    minAlternations    = getOpt(opts, 'MinAlternations', 3);
    chatterMinCoreHrs  = getOpt(opts, 'ChatterMinCoreHours', 1);
    chatterMinFraction = getOpt(opts, 'ChatterMinFraction', 0.20);
    outputDir          = getOpt(opts, 'OutputDir', '');
    verbose            = getOpt(opts, 'Verbose', false);
    sleepOpt           = getOpt(opts, 'Sleep', []);

    if ~isnumeric(dt) || ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error('detect_shuffling:InvalidDt', 'opts.Dt must be a positive scalar.');
    end
    if ~isnumeric(chatterMinFraction) || ~isscalar(chatterMinFraction) || ...
            chatterMinFraction <= 0 || chatterMinFraction > 1
        error('detect_shuffling:InvalidChatterFraction', ...
            'opts.ChatterMinFraction must be in (0, 1].');
    end
    if ~isnumeric(chatterMinCoreHrs) || ~isscalar(chatterMinCoreHrs) || ...
            chatterMinCoreHrs <= 0
        error('detect_shuffling:InvalidChatterHours', ...
            'opts.ChatterMinCoreHours must be a positive scalar.');
    end

    notes = {};

    %% Load and clean

    if isstruct(source) && isfield(source, 'Cleaned')
        pre = source;
        N   = pre.NChannels;
    else
        pre = PreProcessV2(source, N);
    end

    %% Stage 1 QC -- before anything else touches a fly

    qcOpts = struct();
    for f = {'Window', 'DiedMidRecording', 'FlyLabels'}
        if isfield(opts, f{1})
            qcOpts.(f{1}) = opts.(f{1});
        end
    end
    qc = shuffling_qc(pre, qcOpts);

    if ~qc.WindowSupplied
        notes{end+1} = ['No analysis window supplied: whole files used. ', ...
            'Episode indices are NOT aligned to sleep_out.'];
    end

    %% Resolve the sleep trace, if any

    sleep = resolve_sleep(sleepOpt, N, dt);
    if isempty(sleep)
        notes{end+1} = ['No sleep trace supplied: OverlapsSleep and ', ...
            'OverlapsWake are NaN.'];
    elseif ~qc.WindowSupplied
        warning('detect_shuffling:OverlapWithoutWindow', ...
            ['A sleep trace was supplied but no analysis window was. The ', ...
             'overlap columns are being computed against unaligned indices ', ...
             'and are very likely wrong.']);
        notes{end+1} = 'Overlap computed WITHOUT a window: treat as unvalidated.';
    end

    %% Detect, fly by fly

    ep = empty_episodes();
    perFly = struct( ...
        'FlyID',          (1:N)', ...
        'NEpisodes',      zeros(N,1), ...
        'ShuffleSamples', zeros(N,1), ...
        'ShuffleSeconds', zeros(N,1), ...
        'ShuffleFraction', zeros(N,1), ...
        'MaxNRuns',       zeros(N,1), ...
        'MaxCoreSec',     zeros(N,1), ...
        'WindowSeconds',  zeros(N,1), ...
        'SleepSeconds',   nan(N,1), ...
        'WakeSeconds',    nan(N,1));

    for j = 1:N

        w1 = qc.WindowStart(j);
        w2 = qc.WindowEnd(j);
        perFly.WindowSeconds(j) = (w2 - w1 + 1) * dt;

        % Scored sleep and wake time inside this fly's window, needed by
        % the shuffling-time vs wake-time cross-tab in shuffling_stats.
        [perFly.SleepSeconds(j), perFly.WakeSeconds(j)] = ...
            sleep_totals(sleep, j, w1, w2, dt);

        if qc.Excluded(j)
            % Empty / NoValidStart / StuckFromStart: nothing usable.
            perFly.ShuffleFraction(j) = NaN;
            continue
        end

        events = find_shuffles(pre.Cleaned(w1:w2, j), minAlternations, dt);

        if isempty(events)
            continue
        end

        % Coverage mask, not a sum of durations: consecutive episodes can
        % share their boundary run (see find_shuffles), so summing would
        % double-count those samples.
        covered = false(w2 - w1 + 1, 1);

        % Build this fly's block at its known size, then append it once.
        % Appending row by row would recopy the whole table on every
        % episode, and a busy channel can produce tens of thousands.
        nEv = numel(events);
        blk = blank_block(nEv);

        for e = 1:nEv

            evt = events(e);

            % find_shuffles indices are window-relative; convert to raw
            % file rows so every reported index has one unambiguous frame.
            rawStart = evt.startIdx + w1 - 1;
            rawEnd   = evt.endIdx   + w1 - 1;

            covered(evt.startIdx:evt.endIdx) = true;

            [ovSleep, ovWake] = sleep_overlap(sleep, j, rawStart, rawEnd);

            blk.FlyID(e)               = j;
            blk.FlyLabel{e}            = qc.FlyLabel{j};
            blk.StartSample(e)         = rawStart;
            blk.EndSample(e)           = rawEnd;
            blk.StartSampleInWindow(e) = evt.startIdx;
            blk.StartTimeHours(e)      = (evt.startIdx - 1) * dt / 3600;
            blk.DurationSec(e)         = evt.durationSec;
            blk.CoreDurationSec(e)     = evt.coreDurationSec;
            blk.NRuns(e)               = evt.nRuns;
            blk.BeamA(e)               = evt.values(1);
            blk.BeamB(e)               = evt.values(2);
            blk.MaxCoreRunSec(e)       = evt.maxCoreRunSec;
            blk.OverlapsSleep(e)       = ovSleep;
            blk.OverlapsWake(e)        = ovWake;
        end

        ep = append_block(ep, blk);

        perFly.NEpisodes(j)       = numel(events);
        perFly.ShuffleSamples(j)  = sum(covered);
        perFly.ShuffleSeconds(j)  = sum(covered) * dt;
        perFly.ShuffleFraction(j) = sum(covered) / numel(covered);
        perFly.MaxNRuns(j)        = max([events.nRuns]);
        perFly.MaxCoreSec(j)      = max([events.coreDurationSec]);
    end

    %% Stage 2 QC -- ChatteringArtifact
    %
    % The detector's dual role. A single episode of sustained alternation
    % lasting an hour, or covering a fifth of the recording, is a hardware
    % fault rather than behaviour. Judged on core duration so that a fly
    % which simply sat still either side of a brief shuffle is not caught.

    chatterSecThresh = chatterMinCoreHrs * 3600;

    for j = 1:N
        if qc.Excluded(j) || perFly.NEpisodes(j) == 0
            continue
        end
        byTime = perFly.MaxCoreSec(j) >= chatterSecThresh;
        byFrac = perFly.WindowSeconds(j) > 0 && ...
                 (perFly.MaxCoreSec(j) / perFly.WindowSeconds(j)) >= chatterMinFraction;
        if byTime || byFrac
            qc.Status{j}   = 'ChatteringArtifact';
            qc.Excluded(j) = true;
            qc.Detail{j}   = sprintf(['single alternation episode of %.1f min ', ...
                '(%.1f%% of the %.1f h window) -- hardware artifact, not ', ...
                'biological shuffling'], perFly.MaxCoreSec(j)/60, ...
                100*perFly.MaxCoreSec(j)/perFly.WindowSeconds(j), ...
                perFly.WindowSeconds(j)/3600);
        end
    end

    %% Assemble

    result = struct( ...
        'QC',       qc, ...
        'Episodes', ep, ...
        'PerFly',   perFly, ...
        'Params',   struct( ...
            'N',                  N, ...
            'Dt',                 dt, ...
            'MinAlternations',    minAlternations, ...
            'ChatterMinCoreHours', chatterMinCoreHrs, ...
            'ChatterMinFraction', chatterMinFraction, ...
            'WindowSupplied',     qc.WindowSupplied, ...
            'SleepSupplied',      ~isempty(sleep), ...
            'Source',             pre.FilePath), ...
        'Notes',    {notes});

    % Assigned after construction rather than through struct(), which
    % would try to expand a struct-valued field into a struct array.
    result.Sleep = sleep;

    %% Optional outputs

    if ~isempty(outputDir)
        write_tables(result, outputDir, pre.Name);
    end

    if verbose
        print_summary(result);
    end
end

%% ------------------------------------------------------------------ %%

function sleep = resolve_sleep(sleepOpt, N, dt)
% Validate and normalise opts.Sleep into a per-fly cell representation.

    sleep = [];
    if isempty(sleepOpt)
        return
    end
    if ~isstruct(sleepOpt) || ~isscalar(sleepOpt)
        error('detect_shuffling:InvalidSleep', 'opts.Sleep must be a scalar struct.');
    end

    required = {'Trace', 'BinSeconds', 'StartSample', 'Polarity'};
    for i = 1:numel(required)
        if ~isfield(sleepOpt, required{i}) || isempty(sleepOpt.(required{i}))
            error('detect_shuffling:IncompleteSleep', ...
                ['opts.Sleep is missing .%s. All of .Trace, .BinSeconds, ', ...
                 '.StartSample and .Polarity are required.'], required{i});
        end
    end

    % Polarity must be stated: sleeptrace returns 0 = sleep, which is the
    % opposite of the intuitive reading, and guessing inverts everything.
    switch lower(char(sleepOpt.Polarity))
        case 'zeroissleep'
            zeroIsSleep = true;
        case 'oneissleep'
            zeroIsSleep = false;
        otherwise
            error('detect_shuffling:InvalidPolarity', ...
                ['opts.Sleep.Polarity must be ''ZeroIsSleep'' (the ', ...
                 'sleeptrace activity-wave convention) or ''OneIsSleep''.']);
    end

    if iscell(sleepOpt.Trace)
        traces = sleepOpt.Trace(:)';
    elseif isnumeric(sleepOpt.Trace) || islogical(sleepOpt.Trace)
        if size(sleepOpt.Trace, 2) ~= N
            error('detect_shuffling:InvalidSleepTrace', ...
                'opts.Sleep.Trace must have %d columns, one per channel.', N);
        end
        traces = cell(1, N);
        for j = 1:N
            traces{j} = sleepOpt.Trace(:, j);
        end
    else
        error('detect_shuffling:InvalidSleepTrace', ...
            'opts.Sleep.Trace must be numeric, logical, or a cell array.');
    end

    if numel(traces) ~= N
        error('detect_shuffling:InvalidSleepTrace', ...
            'opts.Sleep.Trace must supply %d channels.', N);
    end

    binSeconds  = expandPerFly(sleepOpt.BinSeconds,  N, 'BinSeconds');
    startSample = expandPerFly(sleepOpt.StartSample, N, 'StartSample');

    if any(binSeconds <= 0)
        error('detect_shuffling:InvalidSleepBin', ...
            'opts.Sleep.BinSeconds must be positive.');
    end

    samplesPerBin = binSeconds / dt;
    ragged = find(abs(samplesPerBin - round(samplesPerBin)) > 1e-9, 1);
    if ~isempty(ragged)
        error('detect_shuffling:InvalidSleepBin', ...
            ['opts.Sleep.BinSeconds must be a whole number of samples at ', ...
             'Dt = %g s; got %g s.'], dt, binSeconds(ragged));
    end

    isSleep = cell(1, N);
    for j = 1:N
        v = double(traces{j}(:));
        if zeroIsSleep
            isSleep{j} = (v == 0);
        else
            isSleep{j} = (v ~= 0);
        end
    end

    sleep = struct( ...
        'IsSleep',       {isSleep}, ...
        'SamplesPerBin', round(samplesPerBin), ...
        'StartSample',   startSample);
end

function v = expandPerFly(v, N, name)
    v = double(v(:))';
    if isscalar(v)
        v = repmat(v, 1, N);
    elseif numel(v) ~= N
        error('detect_shuffling:InvalidSleepField', ...
            'opts.Sleep.%s must be scalar or have %d elements.', name, N);
    end
end

function [ovSleep, ovWake] = sleep_overlap(sleep, j, rawStart, rawEnd)
% Score one episode against the fly's sleep/wake trace.
%
% Both raw file row and sleep bin are mapped through the SAME declared
% offset, so the trim that process_LocomToSleep applied is accounted for
% exactly once. An episode spanning a transition returns true for both;
% it is not forced to a single label.

    ovSleep = NaN;
    ovWake  = NaN;

    if isempty(sleep)
        return
    end

    isSleep = sleep.IsSleep{j};
    nBins   = numel(isSleep);
    if nBins == 0
        return
    end

    spb = sleep.SamplesPerBin(j);
    b1  = floor((rawStart - sleep.StartSample(j)) / spb) + 1;
    b2  = floor((rawEnd   - sleep.StartSample(j)) / spb) + 1;

    % Episode entirely outside the scored window: unknown, not false.
    if b2 < 1 || b1 > nBins
        return
    end

    b1 = max(b1, 1);
    b2 = min(b2, nBins);

    span    = isSleep(b1:b2);
    ovSleep = any(span);
    ovWake  = any(~span);
end

function [sleepSec, wakeSec] = sleep_totals(sleep, j, w1, w2, dt)
% Total scored sleep and wake seconds within one fly's analysis window.
%
% Uses the same offset mapping as sleep_overlap, so the cross-tab
% denominator and the per-episode labels cannot disagree about where the
% trim put things.

    sleepSec = NaN;
    wakeSec  = NaN;

    if isempty(sleep)
        return
    end

    isSleep = sleep.IsSleep{j};
    nBins   = numel(isSleep);
    if nBins == 0
        return
    end

    spb = sleep.SamplesPerBin(j);
    b1  = floor((w1 - sleep.StartSample(j)) / spb) + 1;
    b2  = floor((w2 - sleep.StartSample(j)) / spb) + 1;

    if b2 < 1 || b1 > nBins
        return
    end

    b1 = max(b1, 1);
    b2 = min(b2, nBins);

    span      = isSleep(b1:b2);
    secPerBin = spb * dt;
    sleepSec  = sum(span)  * secPerBin;
    wakeSec   = sum(~span) * secPerBin;
end

function blk = blank_block(n)
% One fly's episodes, preallocated to their known count.
    blk = struct( ...
        'FlyID',               zeros(n,1), ...
        'FlyLabel',            {cell(n,1)}, ...
        'StartSample',         zeros(n,1), ...
        'EndSample',           zeros(n,1), ...
        'StartSampleInWindow', zeros(n,1), ...
        'StartTimeHours',      zeros(n,1), ...
        'DurationSec',         zeros(n,1), ...
        'CoreDurationSec',     zeros(n,1), ...
        'NRuns',               zeros(n,1), ...
        'BeamA',               zeros(n,1), ...
        'BeamB',               zeros(n,1), ...
        'MaxCoreRunSec',       zeros(n,1), ...
        'OverlapsSleep',       zeros(n,1), ...
        'OverlapsWake',        zeros(n,1));
end

function ep = append_block(ep, blk)
    names = fieldnames(ep);
    for i = 1:numel(names)
        ep.(names{i}) = [ep.(names{i}); blk.(names{i})];
    end
end

function ep = empty_episodes()
    ep = struct( ...
        'FlyID',               zeros(0,1), ...
        'FlyLabel',            {cell(0,1)}, ...
        'StartSample',         zeros(0,1), ...
        'EndSample',           zeros(0,1), ...
        'StartSampleInWindow', zeros(0,1), ...
        'StartTimeHours',      zeros(0,1), ...
        'DurationSec',         zeros(0,1), ...
        'CoreDurationSec',     zeros(0,1), ...
        'NRuns',               zeros(0,1), ...
        'BeamA',               zeros(0,1), ...
        'BeamB',               zeros(0,1), ...
        'MaxCoreRunSec',       zeros(0,1), ...
        'OverlapsSleep',       zeros(0,1), ...
        'OverlapsWake',        zeros(0,1));
end

function v = getOpt(opts, name, default)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = default;
    end
end

function write_tables(result, outputDir, stem)
% Emit the QC and episode tables as tab-delimited text.

    if exist(outputDir, 'dir') ~= 7
        mkdir(outputDir);
    end
    if isempty(stem)
        stem = 'recording';
    end

    qc = result.QC;
    f  = fopen(fullfile(outputDir, [stem '_shuffling_qc.txt']), 'w');
    fprintf(f, 'FlyID\tFlyLabel\tStatus\tExcluded\tFlags\tStuckValue\tWindowStart\tWindowEnd\tDetail\r\n');
    for i = 1:numel(qc.FlyID)
        fprintf(f, '%d\t%s\t%s\t%d\t%s\t%g\t%d\t%d\t%s\r\n', ...
            qc.FlyID(i), qc.FlyLabel{i}, qc.Status{i}, qc.Excluded(i), ...
            strjoin_compat(qc.Flags{i}, '|'), qc.StuckValue(i), ...
            qc.WindowStart(i), qc.WindowEnd(i), qc.Detail{i});
    end
    fclose(f);

    ep = result.Episodes;
    f  = fopen(fullfile(outputDir, [stem '_shuffling_episodes.txt']), 'w');
    fprintf(f, ['FlyID\tFlyLabel\tStartSample\tEndSample\tStartTimeHours\t', ...
        'DurationSec\tCoreDurationSec\tNRuns\tBeamA\tBeamB\tMaxCoreRunSec\t', ...
        'OverlapsSleep\tOverlapsWake\r\n']);
    for i = 1:numel(ep.FlyID)
        fprintf(f, '%d\t%s\t%d\t%d\t%.4f\t%.2f\t%.2f\t%d\t%g\t%g\t%.2f\t%g\t%g\r\n', ...
            ep.FlyID(i), ep.FlyLabel{i}, ep.StartSample(i), ep.EndSample(i), ...
            ep.StartTimeHours(i), ep.DurationSec(i), ep.CoreDurationSec(i), ...
            ep.NRuns(i), ep.BeamA(i), ep.BeamB(i), ep.MaxCoreRunSec(i), ...
            ep.OverlapsSleep(i), ep.OverlapsWake(i));
    end
    fclose(f);
end

function s = strjoin_compat(c, delim)
    if isempty(c)
        s = '';
        return
    end
    s = c{1};
    for i = 2:numel(c)
        s = [s delim c{i}]; %#ok<AGROW>
    end
end

function print_summary(result)
    qc = result.QC;
    fprintf('\n--- detect_shuffling: %s ---\n', result.Params.Source);
    statuses = unique(qc.Status);
    for i = 1:numel(statuses)
        fprintf('  %-20s %d fly(s)\n', statuses{i}, sum(strcmp(qc.Status, statuses{i})));
    end
    fprintf('  %d episode(s) across %d analysable fly(s)\n', ...
        numel(result.Episodes.FlyID), sum(~qc.Excluded));
    if ~result.Params.WindowSupplied
        fprintf('  NOTE: no analysis window supplied; indices unaligned to sleep_out.\n');
    end
    if ~result.Params.SleepSupplied
        fprintf('  NOTE: no sleep trace supplied; overlap columns are NaN.\n');
    end
    fprintf('\n');
end
