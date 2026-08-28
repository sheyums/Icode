function qc = shuffling_qc(pre, opts)
% SHUFFLING_QC Per-fly status table for the shuffling analysis.
%
%   qc = shuffling_qc(pre)
%   qc = shuffling_qc(pre, opts)
%
%   classifies every channel in one recording before any analysis touches
%   it, so that exclusions are auditable independently of whatever
%   question is being asked that week.
%
%   This is stage 1 of the two-stage QC described in the spec. It covers
%   every status that can be decided from the data alone. The remaining
%   status, ChatteringArtifact, is a property of the detector's own output
%   and is therefore applied afterwards by detect_shuffling, which calls
%   shuffling_qc first and then upgrades the affected rows.
%
%   INPUTS
%   ------
%   pre
%       Output struct from PreProcessV2 (needs .Raw, .Cleaned, .Readable,
%       .NChannels, .NSamples).
%
%   opts
%       Optional struct. Fields, all optional:
%
%         .Window
%             NChannels x 2 matrix of [firstSample lastSample] row indices
%             into the raw file, giving the analysis window for each fly.
%             This is how the entrainment trim (Teliminate) and the death
%             truncation that process_LocomToSleep applies are carried
%             into this tool. Pass [] (default) to analyse whole files,
%             which sets the WindowNotSupplied flag on every fly -- see
%             the caution below.
%
%         .DiedMidRecording
%             NChannels x 1 logical, from process_LocomToSleep's own
%             death detection. Flies marked true are given that status.
%             Default false everywhere. Do not rebuild a parallel dead-fly
%             detector here; that decision belongs upstream.
%
%         .FlyLabels
%             1 x NChannels cell array of names for reporting. Default
%             {'Fly01', ...}.
%
%   OUTPUT
%   ------
%   qc
%       Struct of column vectors, one row per fly:
%           .FlyID        Channel index, 1..NChannels
%           .FlyLabel     Name
%           .Status       One of 'OK', 'Empty', 'NoValidStart',
%                         'StuckFromStart', 'DiedMidRecording'
%                         (detect_shuffling may later set
%                         'ChatteringArtifact')
%           .Excluded     True where nothing usable remains
%           .Flags        Cell of additional non-excluding notes
%           .StuckValue   Beam the channel is stuck on, else NaN
%           .WindowStart  First analysed sample (raw file row)
%           .WindowEnd    Last analysed sample (raw file row)
%           .Detail       Human-readable explanation
%       plus scalar fields .NChannels, .WindowSupplied, .Source.
%
%   STATUS DEFINITIONS
%   ------------------
%   Empty
%       The raw column is exclusively {-1, 0} for the whole file: the
%       channel never reported a position. Excluded. This is a genuinely
%       new check -- the data is finite, valid and present, so neither
%       an empty-file nor a non-finite-data check will catch it.
%
%   NoValidStart
%       The raw column never leaves {0, 1}, so PreProcess's leading-value
%       scan has nothing to seed the forward fill with. Excluded, and
%       handled the same way as Empty. This is the case that makes the
%       original PreProcess throw an index-out-of-bounds; see PreProcessV2.
%
%   StuckFromStart
%       The cleaned position is constant across the entire analysis
%       window. Excluded -- there is no movement to detect. Where two or
%       more flies in the same recording are stuck on the SAME beam, each
%       gets a SharedStuckBeam flag: co-occurrence on one beam is a
%       hardware lead, not a coincidence.
%
%   DiedMidRecording
%       Supplied by the caller from process_LocomToSleep. NOT excluded:
%       the data before truncation is real, and .Window already restricts
%       the analysis to it. shuffling_stats decides whether to pool these
%       flies, and by default it does.
%
%   CAUTION: THE WINDOW IS THE HIGHEST-RISK INPUT
%   ---------------------------------------------
%   process_LocomToSleep trims and truncates before anything else runs,
%   so sample indices in sleep_out do NOT correspond 1:1 to raw file
%   rows. If .Window is omitted, this function analyses whole files and
%   every fly is flagged WindowNotSupplied. That is fine for a standalone
%   look at shuffling, but any comparison against sleep_out made under
%   that flag is silently misaligned by the trim length. Nothing will
%   error. Supply the window.
%
%   See also PREPROCESSV2, DETECT_SHUFFLING, FIND_SHUFFLES.

    %% Validate inputs

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    if ~isstruct(opts) || ~isscalar(opts)
        error('shuffling_qc:InvalidOpts', 'opts must be a scalar struct.');
    end
    requireFields(pre, {'Raw', 'Cleaned', 'Readable', 'NChannels', 'NSamples'}, ...
        'shuffling_qc:InvalidPre', 'pre (expected a PreProcessV2 output struct)');

    N       = pre.NChannels;
    nRows   = pre.NSamples;

    window   = getOpt(opts, 'Window', []);
    died     = getOpt(opts, 'DiedMidRecording', false(N, 1));
    labels   = getOpt(opts, 'FlyLabels', defaultLabels(N));

    windowSupplied = ~isempty(window);

    if windowSupplied
        if ~isnumeric(window) || ~isequal(size(window), [N 2])
            error('shuffling_qc:InvalidWindow', ...
                'opts.Window must be a %d x 2 numeric matrix.', N);
        end
        if any(window(:) ~= fix(window(:))) || any(window(:) < 1) || ...
                any(window(:) > nRows) || any(window(:,2) < window(:,1))
            error('shuffling_qc:InvalidWindow', ...
                ['opts.Window entries must be integer sample indices with ', ...
                 '1 <= start <= end <= %d.'], nRows);
        end
    else
        window = [ones(N,1), repmat(nRows, N, 1)];
    end

    died = logical(died(:));
    if numel(died) ~= N
        error('shuffling_qc:InvalidDied', ...
            'opts.DiedMidRecording must have %d elements.', N);
    end
    if numel(labels) ~= N
        error('shuffling_qc:InvalidLabels', ...
            'opts.FlyLabels must have %d elements.', N);
    end

    %% Classify each channel

    status     = repmat({'OK'}, N, 1);
    detail     = repmat({''},   N, 1);
    flags      = repmat({{}},   N, 1);
    stuckValue = nan(N, 1);
    excluded   = false(N, 1);

    for j = 1:N

        rawCol = pre.Raw(:, j);
        w1     = window(j, 1);
        w2     = window(j, 2);
        cleanW = pre.Cleaned(w1:w2, j);

        if all(rawCol == -1 | rawCol == 0)
            % Never reported a position anywhere in the file.
            status{j}   = 'Empty';
            excluded(j) = true;
            detail{j}   = sprintf(['raw column is exclusively {-1,0} for all ', ...
                '%d samples; channel never reported a position'], nRows);

        elseif ~pre.Readable(j)
            % PreProcess's leading scan cannot terminate on this column.
            status{j}   = 'NoValidStart';
            excluded(j) = true;
            detail{j}   = ['raw column never leaves {0,1}, so the cleaning ', ...
                'forward-fill has no trustworthy value to start from ', ...
                '(this is the case that makes PreProcess throw)'];

        elseif all(cleanW == cleanW(1))
            status{j}     = 'StuckFromStart';
            excluded(j)   = true;
            stuckValue(j) = cleanW(1);
            detail{j}     = sprintf(['cleaned position constant at beam %g ', ...
                'across the whole analysis window (samples %d-%d)'], ...
                cleanW(1), w1, w2);

        elseif died(j)
            status{j} = 'DiedMidRecording';
            detail{j} = sprintf(['death detected upstream; analysed over the ', ...
                'truncated window samples %d-%d only'], w1, w2);

        else
            detail{j} = sprintf('analysed over samples %d-%d', w1, w2);
        end

        if ~windowSupplied
            flags{j}{end+1} = 'WindowNotSupplied';
        end
    end

    %% Shared-beam co-occurrence among StuckFromStart channels
    %
    % Two flies stuck on the same beam in the same recording is a
    % hardware lead worth surfacing, not noise to drop silently.

    isStuck = strcmp(status, 'StuckFromStart');
    stuckIdx = find(isStuck);
    for a = 1:numel(stuckIdx)
        ja    = stuckIdx(a);
        peers = stuckIdx(stuckIdx ~= ja & stuckValue(stuckIdx) == stuckValue(ja));
        if ~isempty(peers)
            flags{ja}{end+1} = 'SharedStuckBeam';
            detail{ja} = sprintf('%s; shares stuck beam %g with channel(s) %s', ...
                detail{ja}, stuckValue(ja), mat2str(peers(:)'));
        end
    end

    %% Assemble

    qc = struct( ...
        'FlyID',          (1:N)', ...
        'FlyLabel',       {labels(:)}, ...
        'Status',         {status}, ...
        'Excluded',       excluded, ...
        'Flags',          {flags}, ...
        'StuckValue',     stuckValue, ...
        'WindowStart',    window(:,1), ...
        'WindowEnd',      window(:,2), ...
        'Detail',         {detail}, ...
        'NChannels',      N, ...
        'WindowSupplied', windowSupplied, ...
        'Source',         pre.FilePath);

    if ~windowSupplied
        warning('shuffling_qc:WindowNotSupplied', ...
            ['No analysis window supplied: using whole files. Episode ', ...
             'sample indices will NOT line up with sleep_out, which is ', ...
             'computed after trimming and truncation. Do not compare ', ...
             'against sleep/wake state under this flag.']);
    end
end

%% ------------------------------------------------------------------ %%

function v = getOpt(opts, name, default)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = default;
    end
end

function requireFields(s, names, errId, what)
    if ~isstruct(s) || ~isscalar(s)
        error(errId, '%s must be a scalar struct.', what);
    end
    for i = 1:numel(names)
        if ~isfield(s, names{i})
            error(errId, '%s is missing required field .%s', what, names{i});
        end
    end
end

function labels = defaultLabels(N)
    labels = cell(1, N);
    for i = 1:N
        labels{i} = sprintf('Fly%02d', i);
    end
end
