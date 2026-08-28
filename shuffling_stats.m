function stats = shuffling_stats(result, opts)
% SHUFFLING_STATS Population-level views of detected shuffling.
%
%   stats = shuffling_stats(result)
%   stats = shuffling_stats(result, opts)
%
%   summarises the output of detect_shuffling across flies. Everything is
%   computed separately at each alternation-length threshold, because a
%   3-run event and a 50-run event are unlikely to be the same phenomenon
%   and pooling them would hide that.
%
%   Re-thresholding here is exactly equivalent to re-running the detector
%   at the higher threshold: find_shuffles extends each alternation as far
%   as it goes and only then applies the minimum, so filtering the
%   recorded episodes on NRuns yields the identical set. Detect once at
%   the lowest threshold of interest and slice here.
%
%   INPUTS
%   ------
%   result
%       Output of detect_shuffling.
%
%   opts
%       Optional struct:
%         .Thresholds      Vector of minimum NRuns values. Default
%                          [3 5 10]. Must be >= the MinAlternations the
%                          detector ran at, or the slice would be
%                          incomplete.
%         .IncludeStatuses Cell of QC statuses to pool. Default
%                          {'OK', 'DiedMidRecording'}. The spec's
%                          population views are defined over Status ==
%                          'OK'; DiedMidRecording is included by default
%                          because that fly's pre-truncation data is
%                          real and .Window already restricts the
%                          analysis to it. Pass {'OK'} for the strict
%                          reading.
%         .TraceBinHours   Bin width for the per-fly time trace meant for
%                          plotting. Default 1 h.
%         .RateBinHours    Bin width for per-fly rates used
%                          statistically. Default 6 h. Deliberately
%                          coarser than the trace: at the incidence rates
%                          seen so far (order 10-30/day at n >= 4), hourly
%                          bins are too sparse per fly to be a stable
%                          rate.
%         .WakeShareFlag   Flag a fly when detected shuffling accounts
%                          for at least this share of its scored wake
%                          time. Default 0.25.
%
%   OUTPUT
%   ------
%   stats
%       Struct with fields:
%         .Thresholds      The thresholds used
%         .ByThreshold     1 x numel(Thresholds) struct array, each with
%                            .Threshold
%                            .FlyID, .FlyLabel
%                            .NEpisodes        per fly
%                            .EpisodesPerDay   per fly
%                            .ShuffleSeconds   per fly, deduplicated
%                            .ShuffleFraction  per fly
%                            .TraceCounts      nFlies x nTraceBins
%                            .TraceEdgesHours
%                            .RateCounts       nFlies x nRateBins
%                            .RatePerHour      nFlies x nRateBins
%                            .RateEdgesHours
%                            .BeamPairs        recurrence table
%                            .WakeCrossTab     see below, or []
%         .IncludedFlyID   Flies pooled
%         .ExcludedSummary Counts by QC status
%         .DayNightIncidence
%                          Not computed. Day/night needs circadian phase
%                          classification, which is out of scope for this
%                          pass; this field carries that statement so it
%                          is visible in the output rather than silently
%                          absent. When it is added it should be paired
%                          per fly (each fly contributing one day rate
%                          and one night rate), not pooled across flies,
%                          since the design is repeated measures.
%         .Notes
%
%   BEAM PAIRS
%   ----------
%   .BeamPairs counts episodes and flies per unordered beam pair. Beam
%   identity is not biologically meaningful, so this is not an analysis
%   axis -- it is a hardware lead. The same pair recurring across
%   independent experiments is worth chasing.
%
%   WAKE CROSS-TAB
%   --------------
%   .WakeCrossTab is populated only when detect_shuffling was given a
%   sleep trace. It is the softer, continuous version of the
%   ChatteringArtifact check: a fly whose shuffling accounts for a large
%   share of its scored wake time is suspect even when no single episode
%   was long enough to trip the hard threshold.
%
%   See also DETECT_SHUFFLING, FIND_SHUFFLES, VALIDATE_SHUFFLING.

    %% Options

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    thresholds      = getOpt(opts, 'Thresholds',      [3 5 10]);
    includeStatuses = getOpt(opts, 'IncludeStatuses', {'OK', 'DiedMidRecording'});
    traceBinHours   = getOpt(opts, 'TraceBinHours',   1);
    rateBinHours    = getOpt(opts, 'RateBinHours',    6);
    wakeShareFlag   = getOpt(opts, 'WakeShareFlag',   0.25);

    if ~isnumeric(thresholds) || isempty(thresholds) || any(thresholds < 2)
        error('shuffling_stats:InvalidThresholds', ...
            'opts.Thresholds must be a non-empty numeric vector of values >= 2.');
    end
    detectedAt = result.Params.MinAlternations;
    if any(thresholds < detectedAt)
        error('shuffling_stats:ThresholdBelowDetection', ...
            ['Threshold %g is below the detector''s MinAlternations (%g). ', ...
             'Episodes shorter than %g were never recorded, so that slice ', ...
             'would be incomplete. Re-run detect_shuffling lower.'], ...
            min(thresholds), detectedAt, detectedAt);
    end
    if traceBinHours <= 0 || rateBinHours <= 0
        error('shuffling_stats:InvalidBin', 'Bin widths must be positive.');
    end
    if rateBinHours < traceBinHours
        warning('shuffling_stats:RateBinTooFine', ...
            ['RateBinHours (%g) is finer than TraceBinHours (%g). Per-fly ', ...
             'rates in bins this narrow are usually too sparse to be stable.'], ...
            rateBinHours, traceBinHours);
    end

    %% Select flies

    qc       = result.QC;
    included = false(numel(qc.FlyID), 1);
    for i = 1:numel(includeStatuses)
        included = included | strcmp(qc.Status, includeStatuses{i});
    end
    flyIdx = find(included);

    if isempty(flyIdx)
        warning('shuffling_stats:NoFliesIncluded', ...
            'No flies match IncludeStatuses; returning empty statistics.');
    end

    windowHours = result.PerFly.WindowSeconds(flyIdx) / 3600;
    maxHours    = max([windowHours(:); traceBinHours; rateBinHours]);

    traceEdges = 0:traceBinHours:ceil(maxHours / traceBinHours) * traceBinHours;
    rateEdges  = 0:rateBinHours:ceil(maxHours / rateBinHours) * rateBinHours;

    ep = result.Episodes;

    %% One slice per alternation-length threshold

    byThreshold = repmat(empty_slice(), 1, numel(thresholds));

    for t = 1:numel(thresholds)

        thr  = thresholds(t);
        keep = ep.NRuns >= thr;

        nFlies      = numel(flyIdx);
        nEpisodes   = zeros(nFlies, 1);
        shuffleSec  = zeros(nFlies, 1);
        traceCounts = zeros(nFlies, max(numel(traceEdges) - 1, 1));
        rateCounts  = zeros(nFlies, max(numel(rateEdges)  - 1, 1));

        for a = 1:nFlies

            j    = flyIdx(a);
            mine = keep & (ep.FlyID == j);

            nEpisodes(a) = sum(mine);

            % Deduplicated: episodes can share a boundary run, so total
            % shuffling time is a coverage measure, not a sum of
            % durations. PerFly.ShuffleSeconds is already deduplicated at
            % the detector's own threshold; above it, rebuild the mask.
            if thr <= detectedAt
                shuffleSec(a) = result.PerFly.ShuffleSeconds(j);
            else
                shuffleSec(a) = covered_seconds( ...
                    ep.StartSample(mine), ep.EndSample(mine), ...
                    result.Params.Dt);
            end

            if any(mine)
                traceCounts(a, :) = bin_counts(ep.StartTimeHours(mine), traceEdges);
                rateCounts(a, :)  = bin_counts(ep.StartTimeHours(mine), rateEdges);
            end
        end

        windowDays = result.PerFly.WindowSeconds(flyIdx) / 86400;
        perDay     = nEpisodes ./ max(windowDays, eps);

        slice = struct( ...
            'Threshold',       thr, ...
            'FlyID',           qc.FlyID(flyIdx), ...
            'FlyLabel',        {qc.FlyLabel(flyIdx)}, ...
            'NEpisodes',       nEpisodes, ...
            'EpisodesPerDay',  perDay, ...
            'ShuffleSeconds',  shuffleSec, ...
            'ShuffleFraction', shuffleSec ./ max(result.PerFly.WindowSeconds(flyIdx), eps), ...
            'TraceCounts',     traceCounts, ...
            'TraceEdgesHours', traceEdges, ...
            'RateCounts',      rateCounts, ...
            'RatePerHour',     rateCounts / rateBinHours, ...
            'RateEdgesHours',  rateEdges, ...
            'BeamPairs',       beam_pairs(ep, keep, flyIdx), ...
            'WakeCrossTab',    wake_cross_tab(result, flyIdx, shuffleSec, wakeShareFlag));

        byThreshold(t) = slice;
    end

    %% Excluded summary

    allStatuses = unique(qc.Status);
    excCount    = zeros(numel(allStatuses), 1);
    for i = 1:numel(allStatuses)
        excCount(i) = sum(strcmp(qc.Status, allStatuses{i}));
    end

    notes = {};
    if ~result.Params.WindowSupplied
        notes{end+1} = ['Windows were not supplied to the detector: elapsed ', ...
            'hours are measured from raw file start, not from the trimmed ', ...
            'recording start.'];
    end
    if ~result.Params.SleepSupplied
        notes{end+1} = 'No sleep trace: WakeCrossTab is empty.';
    end

    stats = struct( ...
        'Thresholds',        thresholds, ...
        'ByThreshold',       byThreshold, ...
        'IncludedFlyID',     qc.FlyID(flyIdx), ...
        'ExcludedSummary',   struct('Status', {allStatuses}, 'Count', excCount), ...
        'DayNightIncidence', ['not computed -- requires circadian phase ', ...
            'classification, deferred from this pass; when added, pair per ', ...
            'fly rather than pooling'], ...
        'Notes',             {notes});
end

%% ------------------------------------------------------------------ %%

function c = bin_counts(values, edges)
% Counts per half-open bin [edges(k), edges(k+1)), last bin closed.

    nBins = max(numel(edges) - 1, 1);
    c     = zeros(1, nBins);
    if isempty(values)
        return
    end
    idx = floor(values / (edges(2) - edges(1))) + 1;
    idx = min(max(idx, 1), nBins);
    for k = 1:numel(idx)
        c(idx(k)) = c(idx(k)) + 1;
    end
end

function sec = covered_seconds(startSamples, endSamples, dt)
% Total distinct samples covered by a set of possibly overlapping spans.
%
% Sorting and merging rather than summing lengths, because episodes may
% share a boundary run.

    sec = 0;
    if isempty(startSamples)
        return
    end
    [s, ord] = sort(startSamples(:));
    e = endSamples(:);
    e = e(ord);

    total    = 0;
    curStart = s(1);
    curEnd   = e(1);
    for k = 2:numel(s)
        if s(k) <= curEnd + 1
            curEnd = max(curEnd, e(k));
        else
            total    = total + (curEnd - curStart + 1);
            curStart = s(k);
            curEnd   = e(k);
        end
    end
    total = total + (curEnd - curStart + 1);
    sec   = total * dt;
end

function bp = beam_pairs(ep, keep, flyIdx)
% Unordered beam-pair recurrence. A hardware lead, not an analysis axis.

    mine = keep;
    if ~isempty(flyIdx)
        inFly = false(size(keep));
        for a = 1:numel(flyIdx)
            inFly = inFly | (ep.FlyID == flyIdx(a));
        end
        mine = keep & inFly;
    end

    lo = min(ep.BeamA(mine), ep.BeamB(mine));
    hi = max(ep.BeamA(mine), ep.BeamB(mine));
    fl = ep.FlyID(mine);

    if isempty(lo)
        bp = struct('BeamLow', zeros(0,1), 'BeamHigh', zeros(0,1), ...
            'NEpisodes', zeros(0,1), 'NFlies', zeros(0,1));
        return
    end

    pairs = unique([lo hi], 'rows');
    nEp   = zeros(size(pairs, 1), 1);
    nFly  = zeros(size(pairs, 1), 1);
    for k = 1:size(pairs, 1)
        m       = (lo == pairs(k,1)) & (hi == pairs(k,2));
        nEp(k)  = sum(m);
        nFly(k) = numel(unique(fl(m)));
    end

    [~, ord] = sort(nEp, 'descend');
    bp = struct( ...
        'BeamLow',   pairs(ord,1), ...
        'BeamHigh',  pairs(ord,2), ...
        'NEpisodes', nEp(ord), ...
        'NFlies',    nFly(ord));
end

function ct = wake_cross_tab(result, flyIdx, shuffleSec, wakeShareFlag)
% Shuffling time as a share of scored wake time, per fly.

    ct = [];
    if ~result.Params.SleepSupplied || isempty(flyIdx)
        return
    end

    wakeSec = result.PerFly.WakeSeconds(flyIdx);
    share   = shuffleSec ./ max(wakeSec, eps);
    share(~isfinite(wakeSec) | wakeSec <= 0) = NaN;

    ct = struct( ...
        'FlyID',              result.PerFly.FlyID(flyIdx), ...
        'ShuffleSeconds',     shuffleSec, ...
        'WakeSeconds',        wakeSec, ...
        'SleepSeconds',       result.PerFly.SleepSeconds(flyIdx), ...
        'ShuffleShareOfWake', share, ...
        'Flagged',            share >= wakeShareFlag, ...
        'FlagThreshold',      wakeShareFlag);
end

function s = empty_slice()
    s = struct( ...
        'Threshold', 0, 'FlyID', [], 'FlyLabel', {{}}, 'NEpisodes', [], ...
        'EpisodesPerDay', [], 'ShuffleSeconds', [], 'ShuffleFraction', [], ...
        'TraceCounts', [], 'TraceEdgesHours', [], 'RateCounts', [], ...
        'RatePerHour', [], 'RateEdgesHours', [], 'BeamPairs', [], ...
        'WakeCrossTab', []);
end

function v = getOpt(opts, name, default)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = default;
    end
end
