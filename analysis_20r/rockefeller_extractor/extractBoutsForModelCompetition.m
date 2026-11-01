function [boutMatrix, meta] = extractBoutsForModelCompetition(sleep_out, phase, boutType, options)
% EXTRACTBOUTSFORMODELCOMPETITION Reshapes pipeline bout output for multiexponential modeling.
%
%   CONTEXT & PURPOSE:
%   'process_LocomToSleep' stores bouts inside a nested '.stats' structure, one
%   fly at a time, with durations in DATA POINTS (bins) and a phase tag on each
%   bout. 'CompareHeterogeneityVsMultiexponential' wants a single NaN-padded
%   [nFlies x maxBouts] numeric matrix. This function does that reshaping and
%   nothing else.
%
%   It does NOT compute xmin. The left-truncation point is a modeling decision
%   that belongs to the caller, and inventing one here would silently couple the
%   fit to assumptions this function has no business making (see XMIN below).
%   The one arithmetic operation applied to the data is the unavoidable scaling
%   from bins to SECONDS (durations are stored in bins; seconds is the unit
%   'process_LocomToSleep' works in). No other unit conversion is offered.
%
%   SYNTAX:
%       [B, meta] = extractBoutsForModelCompetition(S)
%       [B, meta] = extractBoutsForModelCompetition(S, phase)
%       [B, meta] = extractBoutsForModelCompetition(S, phase, boutType)
%       [B, meta] = extractBoutsForModelCompetition(S, phase, boutType, Name=Value)
%
%   EXAMPLE (subjective-night wake bouts):
%       S = process_LocomToSleep(timeInterval=1, segmentMode=3, ...
%                                bridgeWakeSeconds=0, folderPath=D, selectAll=true);
%       [W, meta] = extractBoutsForModelCompetition(S, 'subj_night', 'wake');
%
%       % xmin chosen by YOU, in seconds, from the pipeline's own settings:
%       bridgeBins = floor(S.bridgeWakeSeconds / S.timeInterval);
%       xmin = (bridgeBins + 1) * S.timeInterval;      % = S.timeInterval when unbridged
%
%       % Modeling choices are made HERE, not inside the bridge:
%       R = CompareHeterogeneityVsMultiexponential(W, xmin, 3, ...
%               DistributionType="discrete", ...
%               SamplingInterval=meta.TimeIntervalSeconds);
%
%   INPUTS:
%       sleep_out : (struct) output of process_LocomToSleep. Needs .timeInterval,
%                   .segmentMode, .includedFiles and per-fly .stats.
%       phase     : (char) which bouts to keep, selected by ORIGIN PHASE.
%                   Accepts either the pipeline's field-style names
%                   ('entire','light','dark','subj_day','subj_night') or
%                   segment_sleep_phases' own tag strings
%                   ('Light','Dark','SubjDay','SubjNight').
%                   Default: 'entire' (no filtering).
%       boutType  : (char) 'sleep' or 'wake'. Default: 'sleep'.
%
%   NAME-VALUE:
%       Verbose : (logical) print the extraction report. Default true.
%
%   OUTPUTS:
%       boutMatrix : [nFlies x maxBouts] double, NaN-padded, durations in SECONDS.
%                    Row i corresponds to sleep_out.fly<i> and
%                    sleep_out.includedFiles{i}.
%       meta       : (struct) facts about what was extracted, not derived
%                    modeling quantities:
%                      .Unit               "seconds"
%                      .TimeIntervalSeconds  verbatim sleep_out.timeInterval;
%                                          pass it as the competitor's
%                                          SamplingInterval if you fit discretely
%                      .Phase, .BoutType
%                      .PhaseTag           the tag actually matched; "" when
%                                          phase='entire' (no filtering applied)
%                      .BoutsPerFly, .TotalBouts
%                      .MinObservedSeconds, .MaxObservedSeconds
%                      .OnsetBins          [nFlies x maxBouts] NaN-padded,
%                                          same shape and ordering as boutMatrix:
%                                          OnsetBins(i,j) is the start index of the
%                                          bout whose duration is boutMatrix(i,j),
%                                          verbatim from .stats.*.starts. In BINS,
%                                          not seconds -- it indexes the retained
%                                          trace, so multiplying by TimeIntervalSeconds
%                                          gives elapsed time from the START OF THAT
%                                          TRACE, not from ZT0. NaN is PADDING ONLY:
%                                          a structure carrying durations without
%                                          .starts is rejected, so isnan(OnsetBins) and
%                                          isnan(boutMatrix) agree everywhere.
%                      .DataStartZT        verbatim sleep_out.dataStartZT: the ZT at
%                                          which every .entire trace begins. NaN when
%                                          it is per-fly, in which case the values are
%                                          in sleep_out.fly<i>.dataStartZT.
%                      .OneWholeDay, .PerFlyDayLength   verbatim; day length in bins
%                                          and whether it varies by fly.
%                                          OnsetBins + DataStartZT + OneWholeDay are
%                                          what a day index or a ZT is built from.
%                                          This function does not build either: which
%                                          bin boundary starts 'day 1' is a modeling
%                                          choice, and belongs with xmin at the call site.
%                      .NumCensoredByPipeline  bouts segment_sleep_phases had
%                                          already dropped as boundary-censored.
%                                          A WHOLE-TRACE count: it is not
%                                          restricted to the requested phase,
%                                          because the dropped bouts no longer
%                                          carry a tag to restrict on.
%
%   XMIN -- WHY IT IS NOT RETURNED
%   The competitor left-truncates at xmin and filters with 'row >= xmin', so the
%   value silently determines which bouts enter the likelihood. The correct value
%   depends on facts this function does not own:
%     * sleep: sleeptrace's 5-minute inactivity criterion, i.e.
%           xmin = ceil(300 / S.timeInterval) * S.timeInterval
%       (= 300 exactly whenever timeInterval divides 300).
%     * wake:  ironout absorbs interior wake runs of at most
%           floor(S.bridgeWakeSeconds / S.timeInterval) bins, so
%           xmin = (floor(S.bridgeWakeSeconds / S.timeInterval) + 1) * S.timeInterval.
%       With bridgeWakeSeconds = 0, or below one bin, this is S.timeInterval.
%   Compute xmin as (integer number of bins) * S.timeInterval, the same
%   multiplication used for the durations here. Deriving it any other way can put
%   it one ulp above a minimum-length bout, which deletes every such bout from
%   the fit without warning. meta.MinObservedSeconds lets you check your choice
%   against the data.
%
%   NOTES:
%   1. Durations are read straight from '.stats.*.durations' (bins), not
%      recomputed from starts/ends. Consistency between the two is asserted.
%   2. segment_sleep_phases has ALREADY excluded boundary-censored bouts (those
%      touching the first or last sample). No further censoring is applied here;
%      the count it dropped is reported as meta.NumCensoredByPipeline.
%   3. Bout durations from a binned recording are integer multiples of
%      timeInterval, so what you observe is a bin COUNT, not a duration. The
%      competitor's DistributionType="discrete" branch is the exact grouped-data
%      likelihood for that; its "continuous" default evaluates an exponential
%      density at n*dt and is biased by of order dt/(2*tau) -- negligible for
%      sleep bouts, ~10-20% for few-bin wake bouts. This function does not
%      choose for you -- it reports the bin width and leaves the decision at the
%      call site, with xmin, where it can be stated in a methods section.
%   4. origin_phase is assigned from the bout's STARTING bin and bouts are NOT
%      clipped at phase boundaries -- by design, to avoid fracturing bouts. A
%      bout beginning late in subjective night is therefore returned at full
%      length even if most of it falls in subjective day. This is a property of
%      the data, not of this function, but it acts hardest on the longest bouts,
%      which is where a hyperexponential fit locates its slow component.
%
%   See also PROCESS_LOCOMTOSLEEP, SEGMENT_SLEEP_PHASES,
%            COMPAREHETEROGENEITYVSMULTIEXPONENTIAL.

arguments
    sleep_out struct
    phase (1,:) char {mustBeMember(phase, {'entire', ...
        'light','dark','subj_day','subj_night', ...
        'Light','Dark','SubjDay','SubjNight'})} = 'entire'
    boutType (1,:) char {mustBeMember(boutType, {'sleep','wake'})} = 'sleep'
    options.Verbose (1,1) logical = true
end

%% 0. Structural validation ----------------------------------------------
required = {'timeInterval','includedFiles','segmentMode'};
absent   = required(~isfield(sleep_out, required));
if ~isempty(absent)
    error('extractBouts:NotAPipelineOutput', ...
        ['sleep_out is missing the field(s) %s. This function expects the ' ...
         'structure returned by process_LocomToSleep.'], strjoin(absent, ', '));
end

numFlies = numel(sleep_out.includedFiles);
if numFlies == 0
    error('extractBouts:NoFlies', ...
        'sleep_out.includedFiles is empty -- no flies survived the pipeline.');
end

%% 1. Resolve the requested phase to the tag segment_sleep_phases writes --
% Vocabulary is fixed by segment_sleep_phases' extract_bouts helper:
%   mode 1 -> "Unsegmented" ; mode 2 -> "Light"/"Dark" ; mode 3 -> "SubjDay"/"SubjNight"
tagMap = struct('light',"Light", 'dark',"Dark", ...
                'subj_day',"SubjDay", 'subj_night',"SubjNight", ...
                'Light',"Light", 'Dark',"Dark", ...
                'SubjDay',"SubjDay", 'SubjNight',"SubjNight");

switch sleep_out.segmentMode
    case 1
        allowed  = {'entire'};
        tagsHere = "Unsegmented";
    case 2
        allowed  = {'entire','light','dark','Light','Dark'};
        tagsHere = ["Light" "Dark"];
    case 3
        allowed  = {'entire','subj_day','subj_night','SubjDay','SubjNight'};
        tagsHere = ["SubjDay" "SubjNight"];
    otherwise
        error('extractBouts:UnknownSegmentMode', ...
            'sleep_out.segmentMode = %g is not one of 1, 2 or 3.', sleep_out.segmentMode);
end
if ~ismember(phase, allowed)
    if ismember(phase, {'subj_day','subj_night','SubjDay','SubjNight'})
        needMode = 3;
    else
        needMode = 2;
    end
    error('extractBouts:PhaseNotSegmented', ...
        ['Phase ''%s'' is not available: this structure was built with ' ...
         'segmentMode = %g, under which every bout is tagged %s. Re-run ' ...
         'process_LocomToSleep with segmentMode = %d.'], ...
        phase, sleep_out.segmentMode, strjoin(cellstr(tagsHere), '/'), needMode);
end

filterByPhase = ~strcmp(phase, 'entire');
if filterByPhase
    wantedTag = tagMap.(phase);
else
    wantedTag = "";
end

%% 2. Extract ------------------------------------------------------------
timeIntSec = double(sleep_out.timeInterval);
statsField = sprintf('%s_bouts', boutType);

boutCell  = cell(numFlies, 1);
onsetCell = cell(numFlies, 1);
maxBouts  = 0;
nCensored = 0;
seenTags  = strings(0,1);

for i = 1:numFlies
    flyName = sprintf('fly%d', i);
    if ~isfield(sleep_out, flyName)
        error('extractBouts:MissingFly', ...
            ['includedFiles lists %d flies but field ''%s'' does not exist; ' ...
             'flyN indexing and includedFiles have gone out of step.'], numFlies, flyName);
    end
    flyData = sleep_out.(flyName);

    if ~isfield(flyData, 'stats') || ~isfield(flyData.stats, statsField)
        boutCell{i}  = zeros(1,0);
        onsetCell{i} = zeros(1,0);
        continue
    end
    b = flyData.stats.(statsField);

    durBins = double(b.durations(:));

    % Bout onsets, kept in lockstep with durations through the phase filter.
    % .starts is optional; absent means this fly contributes no day information.
    if ~isfield(b,'starts')
        error('extractBouts:NoOnsets', ...
            ['%s.stats.%s has no ''starts'' field. segment_sleep_phases writes ' ...
             'starts, ends and durations in one block, and duration IS ' ...
             'ends - starts + 1, so a structure with durations but no onsets did ' ...
             'not come from that pipeline.'], flyName, statsField);
    end
    onsetBins = double(b.starts(:));
    if numel(onsetBins) ~= numel(durBins)
        error('extractBouts:OnsetLengthMismatch', ...
            '%s.stats.%s: .starts has %d entries but .durations has %d.', ...
            flyName, statsField, numel(onsetBins), numel(durBins));
    end

    % .durations is authoritative; starts/ends must agree with it. A mismatch
    % means the structure was built by a different version of the segmenter.
    % .starts is guaranteed present by the NoOnsets guard above; .ends is used
    % only for this cross-check, so its absence skips the check rather than
    % rejecting a structure that already carries everything this function needs.
    if isfield(b,'ends')
        recomputed = double(b.ends(:)) - double(b.starts(:)) + 1;
        if numel(recomputed) ~= numel(durBins) || any(recomputed ~= durBins)
            error('extractBouts:InconsistentBoutFields', ...
                ['%s.stats.%s: .durations disagrees with .ends - .starts + 1. ' ...
                 'The bout indexing convention is not what this function assumes.'], ...
                flyName, statsField);
        end
    end

    if isfield(b, 'n_censored')
        nCensored = nCensored + double(b.n_censored);
    end

    if filterByPhase
        if ~isfield(b, 'origin_phase')
            error('extractBouts:NoPhaseTags', ...
                '%s.stats.%s has no ''origin_phase'' field; cannot restrict to ''%s''.', ...
                flyName, statsField, phase);
        end
        tags = string(b.origin_phase(:));
        if numel(tags) ~= numel(durBins)
            error('extractBouts:TagLengthMismatch', ...
                '%s.stats.%s: origin_phase has %d entries but durations has %d.', ...
                flyName, statsField, numel(tags), numel(durBins));
        end
        seenTags  = unique([seenTags; tags]);
        keep      = (tags == wantedTag);
        durBins   = durBins(keep);
        onsetBins = onsetBins(keep);
    end

    boutCell{i}  = durBins(:).' * timeIntSec;     % bins -> seconds
    onsetCell{i} = onsetBins(:).';                % left in BINS, deliberately
    maxBouts = max(maxBouts, numel(boutCell{i}));
end

if maxBouts == 0
    if filterByPhase && ~isempty(seenTags)
        error('extractBouts:PhaseTagNotFound', ...
            ['No %s bout in any fly carries origin_phase = "%s". Tags present ' ...
             'in this structure: %s.'], ...
            boutType, wantedTag, strjoin(cellstr(seenTags(:).'), ', '));
    else
        error('extractBouts:NoBouts', 'No %s bouts were found in any fly.', boutType);
    end
end

%% 3. NaN-padded matrix --------------------------------------------------
boutMatrix = NaN(numFlies, maxBouts);
for i = 1:numFlies
    n = numel(boutCell{i});
    if n > 0
        boutMatrix(i, 1:n) = boutCell{i};
    end
end

onsetMatrix = NaN(numFlies, maxBouts);
for i = 1:numFlies
    n = numel(onsetCell{i});
    if n > 0
        onsetMatrix(i, 1:n) = onsetCell{i};
    end
end

%% 4. Metadata -----------------------------------------------------------
boutsPerFly = cellfun(@numel, boutCell);
allDur      = boutMatrix(isfinite(boutMatrix));

% meta carries two things only: verbatim copies of pipeline metadata, and an
% inventory of what was extracted. No thresholds, no modeling defaults, no
% pre-built argument lists -- those belong to the caller.
meta = struct();
meta.Phase                 = phase;
meta.PhaseTag              = wantedTag;
meta.BoutType              = boutType;
meta.Unit                  = "seconds";
meta.TimeIntervalSeconds   = timeIntSec;          % verbatim sleep_out.timeInterval
meta.BoutsPerFly           = boutsPerFly(:);
meta.TotalBouts            = sum(boutsPerFly);
meta.MinObservedSeconds    = min(allDur);
meta.MaxObservedSeconds    = max(allDur);
meta.NumCensoredByPipeline = nCensored;
meta.OnsetBins             = onsetMatrix;         % bins, same shape as boutMatrix
if isfield(sleep_out,'dataStartZT'),     meta.DataStartZT     = sleep_out.dataStartZT;     else, meta.DataStartZT     = NaN;   end
if isfield(sleep_out,'onewholeday'),     meta.OneWholeDay     = sleep_out.onewholeday;     else, meta.OneWholeDay     = NaN;   end
if isfield(sleep_out,'perFlyDayLength'), meta.PerFlyDayLength = sleep_out.perFlyDayLength; else, meta.PerFlyDayLength = false; end

if options.Verbose
    fprintf('extractBoutsForModelCompetition: %s bouts, phase ''%s'' (tag "%s"), seconds\n', ...
        boutType, phase, wantedTag);
    fprintf('  flies %d | bouts %d (median %g/fly) | %d censored by segmenter\n', ...
        numFlies, meta.TotalBouts, median(boutsPerFly), nCensored);
    fprintf('  observed range %g to %g s (bin = %g s)\n', ...
        meta.MinObservedSeconds, meta.MaxObservedSeconds, timeIntSec);
end
end