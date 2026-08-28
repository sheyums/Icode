function events = find_shuffles(pos, min_alternations, dt)
% FIND_SHUFFLES Detect sustained alternation between two adjacent beams.
%
%   events = find_shuffles(pos)
%   events = find_shuffles(pos, min_alternations)
%   events = find_shuffles(pos, min_alternations, dt)
%
%   finds episodes in which a fly repeatedly alternates between two
%   adjacent beam positions (e.g. 15,14,15,14,15,...), which we call
%   "shuffling". The input is a cleaned position trace at the native
%   acquisition resolution, i.e. the output of PreProcess / PreProcessV2,
%   NOT binned locomotor activity.
%
%   The point of operating here is that sum(abs(diff(position))), as
%   computed by multibeam_locomotion, conflates path length with net
%   displacement by construction. Oscillating in place and walking
%   steadily produce the same activity scalar. The distinction only
%   survives at this resolution.
%
%   INPUTS
%   ------
%   pos
%       Numeric vector of beam-index positions, one sample per row, at the
%       native sampling resolution. Must be finite. Orientation does not
%       matter; it is coerced to a column internally.
%
%   min_alternations
%       Minimum number of runs an episode must contain to be reported.
%       Default 3, i.e. "a genuine excursion and return" (A -> B -> A).
%       Note this counts RUNS, not transitions: a 3-run episode contains
%       2 position changes.
%
%   dt
%       Seconds per sample. Default 0.2 (5 Hz). Used only to convert the
%       sample-count fields into seconds; it does not affect detection.
%
%   OUTPUT
%   ------
%   events
%       1xM struct array, one element per detected episode, with fields:
%
%           .startIdx          Sample index of the first sample of the
%                              episode's first run.
%           .values            Beam values of the episode's runs, in order.
%           .runLengths        Length in samples of each of those runs.
%           .nRuns             Number of runs (numel(.values)).
%           .endIdx            Sample index of the last sample of the
%                              episode's last run.
%           .durationSec       (endIdx - startIdx + 1) * dt. The OUTER
%                              duration -- see the caution below.
%           .coreDurationSec   Duration of the interior runs only, i.e.
%                              the elapsed time between the episode's
%                              first and last position change. Equal to
%                              (runStarts(last) - runStarts(first+1)) * dt.
%           .maxCoreRunSec     Longest interior run, in seconds. Small
%                              values mean dense chatter; large values
%                              mean slow back-and-forth.
%
%       The first four fields are exactly those named in the spec and are
%       unchanged. The remaining four are additive conveniences.
%
%   IMPORTANT: OUTER vs CORE DURATION
%   ---------------------------------
%   The first and last runs of an episode are bounded on their outer side
%   by something that is NOT part of the alternation. A fly that sits on
%   beam 15 for two hours, shuffles 15,14,15, then sits on 15 for another
%   hour produces a single 3-run episode whose first and last runs are
%   enormous. .durationSec would report three hours of "shuffling", which
%   is wrong.
%
%   .coreDurationSec excludes both boundary runs and is therefore the
%   honest measure of how long the fly was actually alternating. Use it
%   for anything quantitative (episode duration statistics, the
%   ChatteringArtifact QC threshold, total shuffling time). .durationSec
%   is retained because .startIdx and .endIdx delimit it, and the
%   sleep/wake overlap check needs that full extent.
%
%   IMPORTANT: EPISODES MAY SHARE A BOUNDARY RUN
%   --------------------------------------------
%   After closing an episode at run j the scan resumes AT j, not at j+1,
%   so a run that terminates one alternation can begin the next. For
%   example [15 16 15 14 15 14] yields a 3-run episode over runs 1-3 and
%   a 4-run episode over runs 3-6, sharing run 3. This is deliberate --
%   the shared run genuinely participates in both patterns -- but it means
%   episode durations MUST NOT simply be summed to get total shuffling
%   time. Build a coverage mask over samples and count that instead;
%   detect_shuffling does this.
%
%   EXAMPLE
%   -------
%       events = find_shuffles([15 16 15 16 15], 3, 0.2);
%       events.nRuns   % 5
%
%   See also DETECT_SHUFFLING, PREPROCESSV2, SHUFFLING_QC.

    %% Defaults

    if nargin < 2 || isempty(min_alternations)
        min_alternations = 3;
    end
    if nargin < 3 || isempty(dt)
        dt = 0.2;
    end

    %% Validate

    if ~isnumeric(pos) && ~islogical(pos)
        error('find_shuffles:InvalidPos', 'pos must be numeric.');
    end
    if ~isvector(pos) && ~isempty(pos)
        error('find_shuffles:InvalidPos', 'pos must be a vector.');
    end
    if ~isnumeric(min_alternations) || ~isscalar(min_alternations) || ...
            min_alternations < 2 || min_alternations ~= fix(min_alternations)
        error('find_shuffles:InvalidMinAlternations', ...
            'min_alternations must be an integer scalar >= 2.');
    end
    if ~isnumeric(dt) || ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error('find_shuffles:InvalidDt', 'dt must be a positive finite scalar.');
    end

    pos = double(pos(:));

    if ~all(isfinite(pos))
        error('find_shuffles:NonFinitePos', ...
            ['pos contains NaN or Inf. Clean the trace first; a NaN would ', ...
             'silently be absorbed into the surrounding run.']);
    end

    events = empty_events();

    if numel(pos) < 2
        return
    end

    %% Decompose the trace into constant-value runs

    changePts  = find(diff(pos) ~= 0) + 1;
    runStarts  = [1; changePts(:)];
    runValues  = pos(runStarts);
    runLengths = diff([runStarts; numel(pos) + 1]);

    R = numel(runValues);

    %% Scan for alternations
    %
    % Grow the output by capacity doubling rather than events(end+1),
    % which would recopy every stored episode on each append. A chattering
    % channel can carry episodes with hundreds of thousands of runs, so
    % those copies are not free. The emitted episodes are identical.

    capacity = 16;
    events   = repmat(event_template(), 1, capacity);
    nEvents  = 0;

    i = 1;
    while i < R
        if abs(runValues(i+1) - runValues(i)) == 1

            A = runValues(i);
            B = runValues(i+1);
            j = i + 1;

            % Extend while runs keep alternating A,B,A,B,...
            % The expected value depends on the parity of the offset from
            % i, not on runValues(i) alone: even offsets are A, odd are B.
            while j + 1 <= R
                if mod(j + 1 - i, 2) == 0
                    expected = A;
                else
                    expected = B;
                end
                if runValues(j+1) == expected
                    j = j + 1;
                else
                    break
                end
            end

            nRuns = j - i + 1;

            if nRuns >= min_alternations
                nEvents = nEvents + 1;
                if nEvents > capacity
                    capacity = capacity * 2;
                    events(capacity) = event_template();
                end
                events(nEvents) = make_event(i, j, nRuns, ...
                    runStarts, runValues, runLengths, dt);
            end

            % Resume AT j, not j+1: see the header note on shared runs.
            i = j;
        else
            i = i + 1;
        end
    end

    if nEvents == 0
        events = empty_events();
    else
        events = events(1:nEvents);
    end
end

%% ------------------------------------------------------------------ %%

function e = make_event(i, j, nRuns, runStarts, runValues, runLengths, dt)
% Assemble one episode record.

    startIdx = runStarts(i);
    endIdx   = runStarts(j) + runLengths(j) - 1;

    % Interior runs i+1 .. j-1. Non-empty whenever nRuns >= 3.
    if nRuns >= 3
        coreDurationSec = (runStarts(j) - runStarts(i+1)) * dt;
        maxCoreRunSec   = max(runLengths(i+1:j-1)) * dt;
    else
        coreDurationSec = 0;
        maxCoreRunSec   = 0;
    end

    e = struct( ...
        'startIdx',        startIdx, ...
        'values',          runValues(i:j), ...
        'runLengths',      runLengths(i:j), ...
        'nRuns',           nRuns, ...
        'endIdx',          endIdx, ...
        'durationSec',     (endIdx - startIdx + 1) * dt, ...
        'coreDurationSec', coreDurationSec, ...
        'maxCoreRunSec',   maxCoreRunSec);
end

function e = event_template()
% Placeholder used only to preallocate; always overwritten before return.

    e = struct( ...
        'startIdx',        0, ...
        'values',          [], ...
        'runLengths',      [], ...
        'nRuns',           0, ...
        'endIdx',          0, ...
        'durationSec',     0, ...
        'coreDurationSec', 0, ...
        'maxCoreRunSec',   0);
end

function e = empty_events()
% 0x0 struct array carrying the full field set.

    e = struct( ...
        'startIdx',        {}, ...
        'values',          {}, ...
        'runLengths',      {}, ...
        'nRuns',           {}, ...
        'endIdx',          {}, ...
        'durationSec',     {}, ...
        'coreDurationSec', {}, ...
        'maxCoreRunSec',   {});
end
