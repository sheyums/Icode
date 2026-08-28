function out = PreProcessV2(source, N)
% PREPROCESSV2 Non-interactive, crash-guarded copy of PreProcess.
%
%   out = PreProcessV2(filepath)
%   out = PreProcessV2(filepath, N)
%   out = PreProcessV2(rawMatrix, N)
%
%   loads and cleans one multibeam raw file, applying exactly the same
%   cleaning rule as PreProcess.m, but callable in a batch loop and
%   without the crash that PreProcess hits on an unusable channel.
%
%   WHY THIS EXISTS
%   ---------------
%   PreProcess.m is treated as read-only production code, so the two
%   changes the shuffling detector needs go in this renamed copy instead
%   of an edit. Both changes are structural, not semantic:
%
%   1. PreProcess selects its input through uigetfile, so it cannot be
%      pointed at a specific file from a script. This version takes a
%      path (or a raw matrix, for testing).
%
%   2. PreProcess's leading-value scan
%
%          while A(index1,j)==0 | A(index1,j)==1
%              index1 = index1+1;
%
%      has no upper bound. A channel whose values never leave {0,1} --
%      exactly the `Empty` case the QC pass in the spec is meant to
%      catch -- walks index1 off the end of the array and throws
%      "index out of bounds" partway through the file, taking the other
%      63 channels down with it. Here the scan is bounded; a channel with
%      no trustworthy leading value is flagged in .Readable and returned
%      uncleaned for shuffling_qc to classify and exclude.
%
%   3. It returns the RAW matrix alongside the cleaned one. The `Empty`
%      check is defined on raw values ({-1,0} only) and the forward-fill
%      destroys that information, so QC cannot run on cleaned data alone.
%
%   The cleaning arithmetic itself is byte-for-byte the same as
%   PreProcess lines 27-38 and is deliberately NOT improved. In
%   particular, both open questions listed as out of scope in the spec
%   are left exactly as they are:
%
%     - the leading run is back-filled with A(index1,j), so a leading -1
%       is treated as a real position rather than a no-read;
%     - a raw 1 is trusted only when the previous CLEANED value is 2
%       (the A(k-1,j)~=2 conditional), which remains unexplained.
%
%   Do not "fix" either one here. If they are ever resolved, that is a
%   separate change with its own validation.
%
%   INPUTS
%   ------
%   source
%       Either a path to a raw tab-delimited .txt file (3 header lines,
%       read with dlmread exactly as PreProcess does), or a numeric raw
%       data matrix that has already been read.
%
%   N
%       Number of fly channels to clean, columns 1..N. Default 64.
%
%   OUTPUT
%   ------
%   out
%       Structure with fields:
%           .Cleaned    Cleaned data matrix. Columns 1..N are cleaned;
%                       any further columns are passed through untouched,
%                       as PreProcess also does.
%           .Raw        The matrix as read, before cleaning.
%           .Readable   1xN logical. False where the channel has no
%                       trustworthy leading value (the bounded-scan case
%                       above). Such columns of .Cleaned equal .Raw.
%           .NChannels  N
%           .NSamples   Number of rows
%           .FilePath   Source path, or '' when a matrix was passed
%           .Name       File name without extension, or ''
%
%   See also PREPROCESS, DETECT_SHUFFLING, SHUFFLING_QC.

    %% Defaults and validation

    if nargin < 2 || isempty(N)
        N = 64;
    end
    if ~isnumeric(N) || ~isscalar(N) || N < 1 || N ~= fix(N)
        error('PreProcessV2:InvalidN', 'N must be a positive integer scalar.');
    end

    filePath = '';
    name     = '';

    if ischar(source) || (exist('isstring', 'builtin') && isstring(source))
        filePath = char(source);
        if exist(filePath, 'file') ~= 2
            error('PreProcessV2:FileNotFound', 'File not found: %s', filePath);
        end
        [~, name] = fileparts(filePath);
        % Same read as PreProcess: tab-delimited, skip 3 header rows.
        A = dlmread(filePath, '\t', 3, 0);
    elseif isnumeric(source)
        A = double(source);
    else
        error('PreProcessV2:InvalidSource', ...
            'source must be a file path or a numeric matrix.');
    end

    if isempty(A)
        error('PreProcessV2:EmptyFile', ...
            'No data rows read from %s.', filePath);
    end
    if size(A, 2) < N
        error('PreProcessV2:TooFewColumns', ...
            'Expected at least %d channel columns, found %d in %s.', ...
            N, size(A, 2), filePath);
    end
    if ~all(all(isfinite(A(:, 1:N))))
        error('PreProcessV2:NonFiniteData', ...
            'Channel columns 1..%d contain NaN or Inf.', N);
    end

    raw    = A;
    nRows  = size(A, 1);
    ok     = true(1, N);

    %% Clean each channel
    %
    % A is mutated in place so that A(k-1,j) below refers to the already
    % CLEANED previous value, not the raw one. PreProcess relies on this
    % and so does the A(k-1,j)~=2 rule; do not hoist a copy out of the
    % loop.

    for j = 1:N

        % Advance past leading samples that carry no trustworthy value.
        % Bounded, unlike the original.
        index1 = 1;
        while index1 <= nRows && (A(index1, j) == 0 || A(index1, j) == 1)
            index1 = index1 + 1;
        end

        if index1 > nRows
            % Whole channel is {0,1}: nothing to seed the fill with. Leave
            % it raw and let QC exclude the fly.
            ok(j) = false;
            continue
        end

        A(1:index1, j) = A(index1, j);

        for k = 2:nRows
            if A(k, j) == 0 || A(k, j) == -1 || ...
                    (A(k, j) == 1 && A(k-1, j) ~= 2)
                A(k, j) = A(k-1, j);
            end
        end
    end

    %% Assemble output

    out = struct( ...
        'Cleaned',   A, ...
        'Raw',       raw, ...
        'Readable',  ok, ...
        'NChannels', N, ...
        'NSamples',  nRows, ...
        'FilePath',  filePath, ...
        'Name',      name);
end
