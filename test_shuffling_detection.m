function test_shuffling_detection()
%% test_shuffling_detection.m
%
% Tests for the shuffling-detection tools: find_shuffles, PreProcessV2,
% shuffling_qc, detect_shuffling, shuffling_stats, validate_shuffling.
%
% Run from the directory containing those files. Each test prints [PASS]
% or [FAIL] with a short description; a summary line is printed at the end.
%
% Tests covered
% -------------
%  1.  Spec sanity check          - [15 16 15 16 15] gives nRuns=5, not 3
%  2.  Long alternation           - a 9-run alternation is not truncated
%  3.  Non-adjacent beams         - 15,17,15,17 is not shuffling
%  4.  min_alternations gate      - 3-run event dropped at a threshold of 4
%  5.  Multi-sample runs          - runs longer than one sample still count
%  6.  Shared boundary run        - overlapping episodes both reported
%  7.  Core vs outer duration     - boundary dwell excluded from the core
%  8.  find_shuffles validation   - bad inputs raise the right identifiers
%  9.  PreProcessV2 equivalence   - cleaning identical to PreProcess
% 10.  PreProcessV2 crash guard   - all-{0,1} column survives, is flagged
% 11.  File round trip            - reads a real 3-header-line .txt
% 12.  QC statuses                - Empty/NoValidStart/StuckFromStart
% 13.  Shared stuck beam          - co-occurrence flagged as hardware lead
% 14.  ChatteringArtifact         - sustained alternation excluded
% 15.  TRIM ALIGNMENT             - elapsed time and sleep bin survive the
%                                   window offset (highest-risk path)
% 16.  Sleep polarity required    - omitting it errors, inverting it flips
% 17.  Sleep/wake spanning        - an episode may be labelled both
% 18.  Threshold equivalence      - filtering == re-detecting higher
% 19.  Threshold below detection  - refused rather than silently partial
% 20.  Coverage deduplication     - shared samples counted once
% 21.  validate_shuffling gaps    - isolated-wake-gap check behaves
% 22.  Unwindowed overlap warns   - misalignment risk is not silent

nPassed = 0;
nFailed = 0;

tests = { ...
    'Spec sanity check (nRuns=5, not 3)',      @t01_spec_sanity; ...
    'Long alternation not truncated',          @t02_long_alternation; ...
    'Non-adjacent beams are not shuffling',    @t03_non_adjacent; ...
    'min_alternations gate',                   @t04_min_gate; ...
    'Runs longer than one sample',             @t05_multi_sample_runs; ...
    'Shared boundary run',                     @t06_shared_run; ...
    'Core vs outer duration',                  @t07_core_duration; ...
    'find_shuffles input validation',          @t08_validation; ...
    'PreProcessV2 matches PreProcess',         @t09_preprocess_equivalence; ...
    'PreProcessV2 crash guard',                @t10_preprocess_guard; ...
    'Raw .txt file round trip',                @t11_file_roundtrip; ...
    'QC statuses',                             @t12_qc_statuses; ...
    'Shared stuck beam flagged',               @t13_shared_stuck; ...
    'ChatteringArtifact excluded',             @t14_chattering; ...
    'TRIM ALIGNMENT survives the offset',      @t15_trim_alignment; ...
    'Sleep polarity must be declared',         @t16_polarity; ...
    'Episode spanning sleep and wake',         @t17_spanning; ...
    'Threshold filtering == re-detecting',     @t18_threshold_equivalence; ...
    'Threshold below detection refused',       @t19_threshold_guard; ...
    'Coverage deduplication',                  @t20_coverage_dedup; ...
    'validate_shuffling wake gaps',            @t21_validate_gaps; ...
    'Unwindowed overlap warns',                @t22_unwindowed_warns};

fprintf('\nRunning tests for shuffling detection...\n\n');

for i = 1:size(tests, 1)
    try
        tests{i, 2}();
        fprintf('[PASS] %2d. %s\n', i, tests{i, 1});
        nPassed = nPassed + 1;
    catch e
        fprintf('[FAIL] %2d. %s: %s\n', i, tests{i, 1}, e.message);
        nFailed = nFailed + 1;
    end
end

fprintf('\n%d passed, %d failed, %d total\n\n', ...
    nPassed, nFailed, nPassed + nFailed);
end

%% ================== find_shuffles ================================== %%

function t01_spec_sanity()
% The bug this replaced always checked whether the next run returned to
% run i's value, which capped every detection at exactly 3 runs.
    e = find_shuffles([15 16 15 16 15], 3, 0.2);
    assertTrue(numel(e) == 1, 'expected exactly one episode');
    assertTrue(e.nRuns == 5, sprintf('nRuns = %d, expected 5', e.nRuns));
    assertTrue(isequal(e.values(:)', [15 16 15 16 15]), 'values wrong');
end

function t02_long_alternation()
    pos = repmat([9 10], 1, 20);          % 40 runs
    e   = find_shuffles(pos, 3, 0.2);
    assertTrue(numel(e) == 1, 'expected one episode');
    assertTrue(e.nRuns == 40, sprintf('nRuns = %d, expected 40', e.nRuns));
end

function t03_non_adjacent()
    e = find_shuffles([15 17 15 17 15], 3, 0.2);
    assertTrue(isempty(e), 'non-adjacent beams must not be reported');
end

function t04_min_gate()
    pos = [5 15 14 15 5];
    assertTrue(numel(find_shuffles(pos, 3, 0.2)) == 1, 'should pass at 3');
    assertTrue(isempty(find_shuffles(pos, 4, 0.2)),    'should fail at 4');
end

function t05_multi_sample_runs()
    pos = [15 15 15 14 14 15 15 14 14 15];   % 5 runs, 2 samples each
    e   = find_shuffles(pos, 3, 0.2);
    assertTrue(numel(e) == 1, 'expected one episode');
    assertTrue(e.nRuns == 5, sprintf('nRuns = %d, expected 5', e.nRuns));
    assertTrue(isequal(e.runLengths(:)', [3 2 2 2 1]), ...
        sprintf('runLengths = %s', mat2str(e.runLengths(:)')));
end

function t06_shared_run()
% [15 16 15 14 15 14]: runs 1-3 alternate 15/16, runs 3-6 alternate 15/14.
% Run 3 belongs to both. Both must be reported.
    e = find_shuffles([15 16 15 14 15 14], 3, 0.2);
    assertTrue(numel(e) == 2, sprintf('expected 2 episodes, got %d', numel(e)));
    assertTrue(e(1).nRuns == 3 && e(2).nRuns == 4, 'run counts wrong');
    assertTrue(e(2).startIdx == 3, 'second episode should start at the shared run');
end

function t07_core_duration()
% A fly that sits still, shuffles briefly, then sits still again must not
% report the sitting as shuffling.
    pos = [repmat(15, 1, 1000), 14, 15, repmat(9, 1, 1000)];
    e   = find_shuffles(pos, 3, 0.2);
    assertTrue(numel(e) == 1, 'expected one episode');
    assertTrue(e.nRuns == 3, 'expected 3 runs');
    assertTrue(abs(e.durationSec - 1002 * 0.2) < 1e-9, ...
        'outer duration should span the whole dwell');
    assertTrue(abs(e.coreDurationSec - 0.2) < 1e-9, ...
        sprintf('core duration = %g, expected 0.2', e.coreDurationSec));
end

function t08_validation()
    assertError(@() find_shuffles([1 NaN 1], 3, 0.2), 'find_shuffles:NonFinitePos');
    assertError(@() find_shuffles([1 2 1], 1, 0.2),   'find_shuffles:InvalidMinAlternations');
    assertError(@() find_shuffles([1 2 1], 3, 0),     'find_shuffles:InvalidDt');
    e = find_shuffles([], 3, 0.2);
    assertTrue(isempty(e), 'empty input should give no episodes');
end

%% ================== PreProcessV2 =================================== %%

function t09_preprocess_equivalence()
% The cleaning must stay byte-identical to PreProcess for every column the
% original can actually process.
    setSeed(11);
    for t = 1:25
        A      = double(randi([-1 20], 300, 6));
        A(1,:) = 5;                       % give the original a valid start
        orig   = originalClean(A, 6);
        v2     = PreProcessV2(A, 6);
        assertTrue(isequal(orig, v2.Cleaned), ...
            sprintf('cleaning diverged on trial %d', t));
        assertTrue(all(v2.Readable), 'all columns should be readable here');
    end
end

function t10_preprocess_guard()
    A       = double(randi([2 20], 200, 3));
    A(:, 2) = 0;
    A(3, 2) = 1;                          % column 2 never leaves {0,1}

    assertError(@() originalClean(A, 3), '');   % original walks off the end

    v2 = PreProcessV2(A, 3);
    assertTrue(isequal(v2.Readable, [true false true]), 'Readable wrong');
    assertTrue(isequal(v2.Cleaned(:,2), A(:,2)), ...
        'unreadable column should be passed through untouched');
end

function t11_file_roundtrip()
    f = [tempname() '.txt'];
    fid = fopen(f, 'w');
    fprintf(fid, 'header line 1\r\nheader line 2\r\nheader line 3\r\n');
    data = [5 5; 15 9; 14 9; 15 9; 14 9; 15 9];
    for r = 1:size(data, 1)
        fprintf(fid, '%d\t%d\r\n', data(r,1), data(r,2));
    end
    fclose(fid);

    cleanup = onCleanupCompat(f);  %#ok<NASGU>

    pre = PreProcessV2(f, 2);
    assertTrue(pre.NSamples == 6, sprintf('read %d rows, expected 6', pre.NSamples));
    assertTrue(isequal(pre.Cleaned(:,1)', data(:,1)'), 'column 1 misread');

    r = quietly(@() detect_shuffling(pre, struct('N', 2)));
    assertTrue(numel(r.Episodes.FlyID) == 1, 'expected one episode');
    assertTrue(r.Episodes.NRuns(1) == 5, ...
        sprintf('NRuns = %d, expected 5', r.Episodes.NRuns(1)));
end

%% ================== QC ============================================= %%

function t12_qc_statuses()
    A = repmat((2:21)', 20, 6);
    A(:, 1) = 0;  A(7, 1) = -1;           % Empty
    A(:, 2) = 0;  A(9, 2) =  1;           % NoValidStart
    A(:, 3) = 12;                         % StuckFromStart

    pre = PreProcessV2(A, 6);
    qc  = quietly(@() shuffling_qc(pre));

    assertTrue(strcmp(qc.Status{1}, 'Empty'),          ['status 1 = ' qc.Status{1}]);
    assertTrue(strcmp(qc.Status{2}, 'NoValidStart'),   ['status 2 = ' qc.Status{2}]);
    assertTrue(strcmp(qc.Status{3}, 'StuckFromStart'), ['status 3 = ' qc.Status{3}]);
    assertTrue(strcmp(qc.Status{4}, 'OK'),             ['status 4 = ' qc.Status{4}]);
    assertTrue(all(qc.Excluded(1:3)), 'first three should be excluded');
    assertTrue(~qc.Excluded(4), 'fourth should not be excluded');
    assertTrue(qc.StuckValue(3) == 12, 'stuck value wrong');
end

function t13_shared_stuck()
    A = repmat((2:21)', 20, 4);
    A(:, 1) = 12;
    A(:, 2) = 12;                         % same beam as fly 1
    A(:, 3) = 17;                         % different beam

    pre = PreProcessV2(A, 4);
    qc  = quietly(@() shuffling_qc(pre));

    assertTrue(hasFlag(qc.Flags{1}, 'SharedStuckBeam'), 'fly 1 not flagged');
    assertTrue(hasFlag(qc.Flags{2}, 'SharedStuckBeam'), 'fly 2 not flagged');
    assertTrue(~hasFlag(qc.Flags{3}, 'SharedStuckBeam'), 'fly 3 wrongly flagged');
end

function t14_chattering()
% Near-continuous alternation is a hardware fault, and the inactivity-based
% death detection upstream structurally cannot see it.
    A       = repmat((2:21)', 500, 3);
    A(:, 2) = repmat([15; 14], 5000, 1);

    r = quietly(@() detect_shuffling(A, struct('N', 3)));

    assertTrue(strcmp(r.QC.Status{2}, 'ChatteringArtifact'), ...
        ['status 2 = ' r.QC.Status{2}]);
    assertTrue(r.QC.Excluded(2), 'chattering fly should be excluded');
    assertTrue(strcmp(r.QC.Status{1}, 'OK'), 'fly 1 should be unaffected');
end

%% ================== Alignment (highest risk) ======================= %%

function t15_trim_alignment()
% The spec's explicit test: a synthetic fly with a known episode at a known
% time, confirmed to survive trimming.
%
% process_LocomToSleep trims and truncates before computing anything, so
% sleep_out indices are not raw file rows. Here the first 1200 samples
% (240 s) are trimmed. The episode sits at raw row 3001, so it is 1800
% samples = 360 s = 0.1 h into the TRIMMED recording. If the offset were
% dropped, elapsed time would come out as 0.1667 h and the sleep bin would
% be 11 rather than 7 -- both wrong, and neither would raise an error.

    dt        = 0.2;
    nSamples  = 6000;
    trimStart = 1201;

    pos = buildBackground(nSamples);
    pos(3001:3005) = [15 14 15 14 15];

    A = [pos(:), buildBackground(nSamples)'];

    % Sleep trace begins at the trimmed start, 60 s bins (300 samples).
    spb    = 300;
    nBins  = floor((nSamples - trimStart + 1) / spb);
    wave   = ones(nBins, 1);              % activity wave: 1 = wake
    wave(7) = 0;                          % bin 7 scored as sleep
    sleepTrace = [wave, ones(nBins, 1)];

    opts = struct( ...
        'N',      2, ...
        'Dt',     dt, ...
        'Window', [trimStart nSamples; trimStart nSamples], ...
        'Sleep',  struct( ...
            'Trace',       sleepTrace, ...
            'BinSeconds',  spb * dt, ...
            'StartSample', trimStart, ...
            'Polarity',    'ZeroIsSleep'));

    r = detect_shuffling(A, opts);

    idx = find(r.Episodes.FlyID == 1);
    assertTrue(numel(idx) == 1, ...
        sprintf('expected exactly one episode on fly 1, got %d', numel(idx)));

    % Raw file row is preserved exactly.
    assertTrue(r.Episodes.StartSample(idx) == 3001, ...
        sprintf('StartSample = %d, expected 3001', r.Episodes.StartSample(idx)));

    % Window-relative index accounts for the trim.
    assertTrue(r.Episodes.StartSampleInWindow(idx) == 3001 - trimStart + 1, ...
        'window-relative index wrong');

    % Elapsed hours are measured from the trimmed start, not the file start.
    expectedHours = (3001 - trimStart) * dt / 3600;
    assertTrue(abs(r.Episodes.StartTimeHours(idx) - expectedHours) < 1e-12, ...
        sprintf('StartTimeHours = %.6f, expected %.6f (untrimmed would be %.6f)', ...
        r.Episodes.StartTimeHours(idx), expectedHours, 3000 * dt / 3600));

    % And the same offset lands the episode in the sleep bin it belongs to.
    assertTrue(r.Episodes.OverlapsSleep(idx) == 1, ...
        'episode should land in the sleep-scored bin 7');
    assertTrue(r.Episodes.OverlapsWake(idx) == 0, ...
        'episode should not touch a wake bin');
end

function t16_polarity()
    [A, opts] = simpleSleepCase();

    bad = opts;
    bad.Sleep = rmfield(bad.Sleep, 'Polarity');
    assertError(@() detect_shuffling(A, bad), 'detect_shuffling:IncompleteSleep');

    bad2 = opts;
    bad2.Sleep.Polarity = 'yes';
    assertError(@() detect_shuffling(A, bad2), 'detect_shuffling:InvalidPolarity');

    % Flipping the declared convention must flip every label.
    r1 = detect_shuffling(A, opts);
    o2 = opts; o2.Sleep.Polarity = 'OneIsSleep';
    r2 = detect_shuffling(A, o2);

    assertTrue(r1.Episodes.OverlapsSleep(1) ~= r2.Episodes.OverlapsSleep(1), ...
        'polarity change did not flip the sleep label');
end

function t17_spanning()
% An episode straddling a transition touches both states and must not be
% forced to a single label.
    dt   = 0.2;
    n    = 1200;
    pos  = buildBackground(n);
    pos(598:603) = [15 14 15 14 15 14];   % straddles the 600-sample boundary

    spb   = 300;
    wave  = [1; 0; 1; 1];                 % bins: wake, sleep, wake, wake
    opts  = struct('N', 1, 'Dt', dt, ...
        'Window', [1 n], ...
        'Sleep', struct('Trace', wave, 'BinSeconds', spb * dt, ...
            'StartSample', 1, 'Polarity', 'ZeroIsSleep'));

    r = detect_shuffling(pos(:), opts);
    i = find(r.Episodes.FlyID == 1, 1);
    assertTrue(~isempty(i), 'no episode detected');
    assertTrue(r.Episodes.OverlapsSleep(i) == 1 && r.Episodes.OverlapsWake(i) == 1, ...
        sprintf('expected both true, got sleep=%g wake=%g', ...
        r.Episodes.OverlapsSleep(i), r.Episodes.OverlapsWake(i)));
end

%% ================== Stats ========================================== %%

function t18_threshold_equivalence()
% Filtering recorded episodes on NRuns must equal detecting at that
% threshold, otherwise the per-threshold views in the stats are wrong.
    setSeed(3);
    pos = buildBackground(4000);
    pos(101:106)   = [15 14 15 14 15 14];   % 6 runs
    pos(1501:1503) = [8 9 8];               % 3 runs
    pos(2501:2510) = repmat([4 5], 1, 5);   % 10 runs

    at3 = find_shuffles(pos, 3, 0.2);
    at5 = find_shuffles(pos, 5, 0.2);

    filtered = at3([at3.nRuns] >= 5);
    assertTrue(numel(filtered) == numel(at5), ...
        sprintf('filtered %d vs detected %d', numel(filtered), numel(at5)));
    assertTrue(isequal([filtered.startIdx], [at5.startIdx]), ...
        'filtered and re-detected episodes differ');
end

function t19_threshold_guard()
    A = repmat((2:21)', 50, 2);
    A(101:106, 1) = [15 14 15 14 15 14];
    r = quietly(@() detect_shuffling(A, struct('N', 2, 'MinAlternations', 5)));
    assertError(@() shuffling_stats(r, struct('Thresholds', 3)), ...
        'shuffling_stats:ThresholdBelowDetection');
end

function t20_coverage_dedup()
% Episodes sharing a boundary run must not have their overlap counted
% twice in total shuffling time.
    n   = 2000;
    pos = buildBackground(n);
    pos(501:506) = [15 16 15 14 15 14];   % two overlapping episodes

    r = quietly(@() detect_shuffling(pos(:), struct('N', 1, 'Window', [1 n])));

    assertTrue(r.PerFly.NEpisodes(1) == 2, ...
        sprintf('expected 2 episodes, got %d', r.PerFly.NEpisodes(1)));

    summed = sum(r.Episodes.DurationSec);
    assertTrue(r.PerFly.ShuffleSeconds(1) < summed, ...
        'coverage should be strictly less than the sum of durations');
    assertTrue(abs(r.PerFly.ShuffleSeconds(1) - 6 * 0.2) < 1e-9, ...
        sprintf('coverage = %g s, expected %g s', ...
        r.PerFly.ShuffleSeconds(1), 6 * 0.2));
end

%% ================== Validation ===================================== %%

function t21_validate_gaps()
    dt  = 0.2;
    spb = 300;

    % 10 sleep bins, one wake bin, 10 sleep bins: a single isolated gap
    % of exactly one bin, flanked by runs well over 5 minutes.
    wave = [zeros(10,1); 1; zeros(10,1)];

    n    = spb * numel(wave);
    pos  = buildBackground(n);
    pos(501:505) = [15 14 15 14 15];

    opts = struct('N', 1, 'Dt', dt, 'Window', [1 n], ...
        'Sleep', struct('Trace', wave, 'BinSeconds', spb * dt, ...
            'StartSample', 1, 'Polarity', 'ZeroIsSleep'));

    r   = detect_shuffling(pos(:), opts);
    rep = validate_shuffling(r);

    assertTrue(rep.WakeGaps.Ran, 'wake-gap check should have run');
    assertTrue(rep.WakeGaps.NGaps == 1, ...
        sprintf('expected 1 qualifying gap, got %d', rep.WakeGaps.NGaps));
    assertTrue(rep.WakeGaps.NGapsLongerThanOneSample == 0, 'gap should be one bin');
    assertTrue(rep.WakeGaps.Passed, 'check should pass');

    % A 3-run shuffle at 5 Hz spans more than one second, so bridging at
    % the default 1 s cannot absorb it.
    assertTrue(rep.Bridging.Ran && rep.Bridging.Passed, ...
        sprintf('bridging check: %s', rep.Bridging.Reason));
end

function t22_unwindowed_warns()
% Comparing against sleep_out without the trim offset is the silent
% failure the spec warns about, so it must at least be loud.
    [A, opts] = simpleSleepCase();
    opts = rmfield(opts, 'Window');

    lastwarn('');
    warnState = warning('off', 'all');
    r = detect_shuffling(A, opts);
    warning(warnState);

    assertTrue(~r.Params.WindowSupplied, 'window should be absent');
    assertTrue(any(cellfun(@(s) ~isempty(strfind(s, 'WITHOUT a window')), r.Notes)), ...
        'notes should record the unaligned overlap');
end

%% ================== Helpers ======================================== %%

function pos = buildBackground(n)
% Background that varies (so it is not StuckFromStart) but whose values are
% never adjacent to each other or to the injected beams 14/15/16.
    pos = zeros(1, n);
    for k = 1:n
        if mod(floor((k - 1) / 500), 2) == 0
            pos(k) = 5;
        else
            pos(k) = 9;
        end
    end
end

function [A, opts] = simpleSleepCase()
    dt  = 0.2;
    n   = 1200;
    pos = buildBackground(n);
    pos(301:305) = [15 14 15 14 15];

    spb  = 300;
    wave = [1; 0; 1; 1];       % bin 2 (samples 301-600) scored as sleep

    A    = pos(:);
    opts = struct('N', 1, 'Dt', dt, 'Window', [1 n], ...
        'Sleep', struct('Trace', wave, 'BinSeconds', spb * dt, ...
            'StartSample', 1, 'Polarity', 'ZeroIsSleep'));
end

function A = originalClean(A, N)
% PreProcess.m lines 27-38, verbatim, for the equivalence test.
    for j = 1:N
        index1 = 1;
        while A(index1,j)==0 || A(index1,j)==1
            index1 = index1+1;
        end
        A(1:index1,j) = A(index1,j);
        for k = 2:size(A,1)
            if A(k,j)==0 || A(k,j)==-1 || (A(k,j)==1 && A(k-1,j)~=2)
                A(k,j) = A(k-1,j);
            end
        end
    end
end

function out = quietly(fn)
    warnState = warning('off', 'all');
    try
        out = fn();
    catch e
        warning(warnState);
        rethrow(e);
    end
    warning(warnState);
end

function tf = hasFlag(flags, name)
    tf = false;
    for i = 1:numel(flags)
        if strcmp(flags{i}, name)
            tf = true;
            return
        end
    end
end

function assertTrue(cond, msg)
    if ~cond
        error('test:Failed', '%s', msg);
    end
end

function assertError(fn, expectedId)
% Asserts fn errors. When expectedId is non-empty, the identifier must match.
    threw = false;
    try
        warnState = warning('off', 'all');
        fn();
        warning(warnState);
    catch e
        warning(warnState);
        threw = true;
        if ~isempty(expectedId) && ~strcmp(e.identifier, expectedId)
            error('test:Failed', 'expected %s, got %s (%s)', ...
                expectedId, e.identifier, e.message);
        end
    end
    if ~threw
        error('test:Failed', 'expected an error, none raised');
    end
end

function setSeed(s)
    if exist('rng', 'builtin') || exist('rng', 'file')
        rng(s);
    else
        rand('seed', s); randn('seed', s);
    end
end

function c = onCleanupCompat(f)
    if exist('onCleanup', 'class') || exist('onCleanup', 'file')
        c = onCleanup(@() delete(f));
    else
        c = [];
    end
end
