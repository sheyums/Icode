function report = validate_shuffling(result, opts)
% VALIDATE_SHUFFLING Post-build checks on a shuffling detection run.
%
%   report = validate_shuffling(result)
%   report = validate_shuffling(result, opts)
%
%   runs the two checks the spec asks for once the detector exists. Both
%   need the sleep/wake trace; each is skipped, with a reason recorded,
%   when its inputs are absent. Nothing here modifies result.
%
%   CHECK 1: ISOLATED WAKE GAPS
%   ---------------------------
%   Re-runs the earlier finding that a wake gap sandwiched between two
%   sleep runs, each independently long enough to qualify as sleep
%   (>= 5 min by default), is always exactly one sample long. Restricted
%   to the analysable fly subset, since part of the original result may
%   have come from channels that have since been flagged. If the finding
%   no longer holds, that is itself informative -- it means some of the
%   original events came from a channel this QC pass now excludes.
%
%   CHECK 2: BRIDGING CANNOT SWALLOW REAL SHUFFLING
%   -----------------------------------------------
%   Confirms that genuine multi-run shuffling is essentially never
%   affected by ironout's wake-bridging. ironout merges wake gaps shorter
%   than bridgeWakeSeconds (default 1 s) into the surrounding sleep, so an
%   episode is only at risk if it falls entirely inside a wake run that
%   short. This is expected to be empty given the mechanism -- a 3-run
%   alternation at 5 Hz already spans more than a second -- but the spec
%   asks for it to be asserted rather than assumed, so it is measured
%   directly rather than argued.
%
%   Note this reasons about what bridging WOULD do from the wake-run
%   lengths and the threshold; it does not call ironout. That keeps the
%   check runnable without reaching into the sleep pipeline, at the cost
%   of depending on bridgeWakeSeconds being passed correctly.
%
%   INPUTS
%   ------
%   result
%       Output of detect_shuffling. Uses .Sleep, .Episodes, .QC, .PerFly.
%
%   opts
%       Optional struct:
%         .MinSleepMinutes     Sleep-run length that qualifies as sleep
%                              on each side of a gap. Default 5.
%         .BridgeWakeSeconds   ironout's bridging threshold. Default 1.
%         .MinRuns             Episode threshold for check 2. Default 3.
%         .IncludeStatuses     QC statuses to consider. Default
%                              {'OK', 'DiedMidRecording'}.
%
%   OUTPUT
%   ------
%   report
%       Struct with fields .WakeGaps and .Bridging, each carrying .Ran,
%       .Passed, .Skipped, .Reason and check-specific detail, plus
%       .Summary, a printable multi-line string.
%
%   See also DETECT_SHUFFLING, SHUFFLING_STATS.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    minSleepMinutes   = getOpt(opts, 'MinSleepMinutes',   5);
    bridgeWakeSeconds = getOpt(opts, 'BridgeWakeSeconds', 1);
    minRuns           = getOpt(opts, 'MinRuns',           3);
    includeStatuses   = getOpt(opts, 'IncludeStatuses',   {'OK', 'DiedMidRecording'});

    qc       = result.QC;
    included = false(numel(qc.FlyID), 1);
    for i = 1:numel(includeStatuses)
        included = included | strcmp(qc.Status, includeStatuses{i});
    end
    flyIdx = find(included);

    hasSleep = isfield(result, 'Sleep') && ~isempty(result.Sleep);

    report = struct( ...
        'WakeGaps', check_wake_gaps(result, flyIdx, hasSleep, minSleepMinutes), ...
        'Bridging', check_bridging(result, flyIdx, hasSleep, bridgeWakeSeconds, minRuns));

    report.Summary = format_summary(report);
end

%% ------------------------------------------------------------------ %%

function c = check_wake_gaps(result, flyIdx, hasSleep, minSleepMinutes)

    c = struct('Ran', false, 'Passed', false, 'Skipped', true, ...
        'Reason', '', 'NGaps', 0, 'NGapsLongerThanOneSample', 0, ...
        'GapLengthBins', zeros(0,1), 'FlyID', zeros(0,1), ...
        'MinSleepMinutes', minSleepMinutes);

    if ~hasSleep
        c.Reason = 'no sleep trace supplied to detect_shuffling';
        return
    end
    if isempty(flyIdx)
        c.Reason = 'no analysable flies';
        return
    end

    dt        = result.Params.Dt;
    gapLens   = [];
    gapFly    = [];

    for a = 1:numel(flyIdx)

        j       = flyIdx(a);
        isSleep = result.Sleep.IsSleep{j};
        if numel(isSleep) < 3
            continue
        end

        secPerBin  = result.Sleep.SamplesPerBin(j) * dt;
        minSleepBins = (minSleepMinutes * 60) / secPerBin;

        [runVal, runLen] = run_length_encode(isSleep(:));

        % A wake run qualifies when both neighbours are sleep runs that
        % are each independently long enough to count as sleep.
        for k = 2:numel(runVal) - 1
            if ~runVal(k) && runVal(k-1) && runVal(k+1) && ...
                    runLen(k-1) >= minSleepBins && runLen(k+1) >= minSleepBins
                gapLens(end+1,1) = runLen(k);  %#ok<AGROW>
                gapFly(end+1,1)  = j;          %#ok<AGROW>
            end
        end
    end

    c.Ran                      = true;
    c.Skipped                  = false;
    c.NGaps                    = numel(gapLens);
    c.GapLengthBins            = gapLens;
    c.FlyID                    = gapFly;
    c.NGapsLongerThanOneSample = sum(gapLens > 1);
    c.Passed                   = (c.NGaps > 0) && (c.NGapsLongerThanOneSample == 0);

    if c.NGaps == 0
        c.Reason = 'no qualifying wake gaps found in the analysable subset';
    elseif ~c.Passed
        c.Reason = sprintf(['%d of %d qualifying gaps are longer than one ', ...
            'sample; the earlier finding does not hold on this subset'], ...
            c.NGapsLongerThanOneSample, c.NGaps);
    end
end

function c = check_bridging(result, flyIdx, hasSleep, bridgeWakeSeconds, minRuns)

    c = struct('Ran', false, 'Passed', false, 'Skipped', true, ...
        'Reason', '', 'NEpisodesChecked', 0, 'NAtRisk', 0, ...
        'AtRiskEpisodeIdx', zeros(0,1), ...
        'BridgeWakeSeconds', bridgeWakeSeconds, 'MinRuns', minRuns);

    if ~hasSleep
        c.Reason = 'no sleep trace supplied to detect_shuffling';
        return
    end

    ep = result.Episodes;
    if isempty(ep.FlyID)
        c.Ran = true; c.Skipped = false; c.Passed = true;
        c.Reason = 'no episodes to check';
        return
    end

    dt     = result.Params.Dt;
    atRisk = zeros(0,1);
    nCheck = 0;

    for e = 1:numel(ep.FlyID)

        j = ep.FlyID(e);
        if ~any(flyIdx == j) || ep.NRuns(e) < minRuns
            continue
        end

        isSleep = result.Sleep.IsSleep{j};
        if isempty(isSleep)
            continue
        end

        nCheck    = nCheck + 1;
        spb       = result.Sleep.SamplesPerBin(j);
        secPerBin = spb * dt;

        b1 = floor((ep.StartSample(e) - result.Sleep.StartSample(j)) / spb) + 1;
        b2 = floor((ep.EndSample(e)   - result.Sleep.StartSample(j)) / spb) + 1;
        b1 = max(b1, 1);
        b2 = min(b2, numel(isSleep));
        if b2 < b1
            continue
        end

        % At risk only if the whole episode sits inside a single wake run
        % short enough for bridging to absorb it.
        if all(~isSleep(b1:b2))
            [runVal, runLen, runStart] = run_length_encode(isSleep(:));
            k = find(runStart <= b1, 1, 'last');
            if ~isempty(k) && ~runVal(k) && ...
                    (runStart(k) + runLen(k) - 1) >= b2 && ...
                    (runLen(k) * secPerBin) <= bridgeWakeSeconds
                atRisk(end+1,1) = e; %#ok<AGROW>
            end
        end
    end

    c.Ran              = true;
    c.Skipped          = false;
    c.NEpisodesChecked = nCheck;
    c.NAtRisk          = numel(atRisk);
    c.AtRiskEpisodeIdx = atRisk;
    c.Passed           = (c.NAtRisk == 0);

    if ~c.Passed
        c.Reason = sprintf(['%d of %d episodes lie entirely inside a wake ', ...
            'run of <= %g s and could be absorbed by bridging'], ...
            c.NAtRisk, nCheck, bridgeWakeSeconds);
    end
end

function [vals, lens, starts] = run_length_encode(x)
    x      = x(:);
    starts = [1; find(diff(x) ~= 0) + 1];
    vals   = x(starts);
    lens   = diff([starts; numel(x) + 1]);
end

function s = format_summary(report)
    lines = {};
    lines{end+1} = 'validate_shuffling';
    lines{end+1} = sprintf('  wake gaps : %s', verdict(report.WakeGaps));
    if report.WakeGaps.Ran
        lines{end+1} = sprintf('              %d qualifying gaps, %d longer than one sample', ...
            report.WakeGaps.NGaps, report.WakeGaps.NGapsLongerThanOneSample);
    end
    lines{end+1} = sprintf('  bridging  : %s', verdict(report.Bridging));
    if report.Bridging.Ran
        lines{end+1} = sprintf('              %d episodes checked, %d at risk', ...
            report.Bridging.NEpisodesChecked, report.Bridging.NAtRisk);
    end
    s = strjoin_lines(lines);
end

function v = verdict(c)
    if c.Skipped
        v = sprintf('SKIPPED (%s)', c.Reason);
    elseif c.Passed
        v = 'PASS';
    else
        v = sprintf('FAIL (%s)', c.Reason);
    end
end

function s = strjoin_lines(c)
    s = '';
    for i = 1:numel(c)
        s = [s c{i} sprintf('\n')]; %#ok<AGROW>
    end
end

function v = getOpt(opts, name, default)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = default;
    end
end
