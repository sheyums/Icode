function O = o1a_run(csvFile, outDir)
%O1A_RUN  PRESPEC O1-A (beam_hmm/PRESPEC.md, 2026-11-19): sleep- and wake-bout hazards and model fits for the CS LD
%pilot (weeks 2-3, primary set), with the repo-root bout library (BoutHazard, CompareBoutModels).
%
%   Run from the Icode REPO ROOT in MATLAB R2025a+ (not Octave; CLAUDE.md):
%     addpath('beam_hmm');
%     O = o1a_run('beam_hmm/results/o1a/bouts_weeks2-3.csv', 'beam_hmm/results/o1a');
%   Input: the CSV written by beam_hmm/o1.py (columns ch, week, kind, start_s, length_s, start_zt).
%   Output: O (struct), outDir/o1a_results.mat, outDir/o1a_report.txt (the full Command Window output), and
%   hazard figures outDir/o1a_hazard_<kind>_<nb>bins.png.
%   Order (PRESPEC): CompareBoutModels first only to get the best hyperexponential as the hazard NULL; the hazard
%   test is reported first; then the model table. Secondary: hazards by start phase (light ZT 0-12 / dark),
%   DESCRIPTIVE, with the censoring caveat (bouts assigned by start, not clipped).
%   Run time: CompareBoutModels on a few thousand bouts can take tens of minutes per bout type.

assert(isequal(fileparts(which('CompareBoutModels')), pwd), ...
       'Run from the Icode repo root, so the repo copies of the bout functions are used (CLAUDE.md).');
if ~exist(outDir, 'dir'), mkdir(outDir); end
logFile = fullfile(outDir, 'o1a_report.txt'); if exist(logFile, 'file'), delete(logFile); end
diary(logFile); cleanup = onCleanup(@() diary('off'));
T = readtable(csvFile, 'TextType', 'string');
fprintf('O1-A  %s  (%s)\n', csvFile, char(datetime('now')));
fprintf('flies: %s\n', mat2str(unique(T.ch)'));
seed = 20261119; NB = [16 20];
spec = struct('kind', {"wake", "sleep"}, 'xmin', {2, 300});           % xmin: CLAUDE.md formulas (bridging 1 s; 300 s)
O = struct();
for s = 1:numel(spec)
    k = spec(s).kind; xmin = spec(s).xmin;
    sel = T.kind == k; x = double(T.length_s(sel)); zt = T.start_zt(sel);
    fprintf('\n================ %s bouts: n = %d, xmin = %g s, dt = 1 s, median %g s, max %g s\n', ...
            upper(k), numel(x), xmin, median(x), max(x));
    assert(all(x >= xmin), 'a bout below xmin: extraction and xmin disagree');
    % ---- null handle (step 2 run first; its table is printed after the hazard, as pre-specified)
    R = CompareBoutModels(x, xmin, SamplingInterval=1, MaxComponents=5, RandomSeed=seed, Plot=false, Verbose=false);
    tab = R.Table; mdl = string(tab.Model); deg = logical(tab.Degenerate);
    iH = find(startsWith(mdl, "hyperexp") & ~deg, 1);
    assert(~isempty(iH), 'no identified hyperexponential row: no null for the hazard test');
    nullS = R.Fits(iH).SurvivalHandle;
    fprintf('hazard null: %s (best identified hyperexponential by AICc)\n', mdl(iH));
    % ---- step 1: hazard test at the pre-specified bin counts
    for nb = NB
        H = BoutHazard(x, xmin, SamplingInterval=1, NumBins=nb, MinAtRisk=10, NullSurvivalHandle=nullS, ...
                       NullReplicates=999, Bootstrap=999, RandomSeed=seed, SurvivalHandle=nullS, Plot=true);
        fprintf('\n-- %s hazard, %d bins: RiseNullP = %.4f (NullDropped %d), RiseRatio %.3f, trough at %.4g s, peak at %.4g s\n', ...
                k, nb, H.RiseNullP, H.NullDropped, H.RiseRatio, H.BinCenters(H.TroughIndex), H.BinCenters(H.PeakIndex));
        fprintf('   %10s %8s %7s %10s %10s %10s %10s\n', 'center_s', 'atrisk', 'events', 'hazard', 'lower', 'upper', 'null');
        for j = 1:numel(H.Hazard)
            fprintf('   %10.4g %8d %7d %10.3e %10.3e %10.3e %10.3e\n', H.BinCenters(j), H.AtRisk(j), H.Events(j), ...
                    H.Hazard(j), H.Lower(j), H.Upper(j), H.ModelHazard(j));
        end
        fh = gcf; try theme(fh, 'light'); catch, end
        print(fh, fullfile(outDir, sprintf('o1a_hazard_%s_%dbins.png', k, nb)), '-dpng', '-r150'); close(fh);
        O.(k).(sprintf('H%d', nb)) = H;
    end
    p = [O.(k).H16.RiseNullP, O.(k).H20.RiseNullP]; okNull = O.(k).H16.NullDropped == 0 && O.(k).H20.NullDropped == 0;
    if ~okNull, v = 'NOT EVALUABLE (NullDropped > 0)';
    elseif all(p < 0.05), v = 'SUPPORTED (RiseNullP < 0.05 at 16 AND 20 bins)';
    elseif all(p >= 0.05), v = 'NOT SUPPORTED (RiseNullP >= 0.05 at both)';
    else, v = 'MIXED (no claim)'; end
    fprintf('\n>> %s hazard rise: %s\n', upper(k), v); O.(k).verdict = v;
    % ---- step 2 report
    fprintf('\n-- %s model comparison (read in this order):\n', k);
    fprintf('R.Skipped:\n'); disp(R.Skipped);
    disp(tab(:, intersect({'Model', 'k', 'LogLik', 'AICc', 'dAICc', 'AkaikeWeight', 'BIC', 'dBIC', 'Degenerate', 'Reason', 'ParamText'}, ...
                          tab.Properties.VariableNames, 'stable')));
    if isempty(R.Selected), selName = "(none passed the fit test)"; else, selName = string(R.Selected.Model); end
    fprintf('R.Selected: %s\n', selName); disp(R.SelectionPath);   % R.Selected is a table-row struct, not a string
    if ~isempty(R.SelectedFit) && isfield(R.SelectedFit, 'Full') && isfield(R.SelectedFit.Full, 'Diagnostics')
        disp(R.SelectedFit.Full.Diagnostics);
    end
    if isfield(R, 'GoF'), disp(R.GoF); end
    if isfield(R, 'OrderLRT') && ~isempty(R.OrderLRT), fprintf('OrderLRT rungs: %d\n', numel(R.OrderLRT)); end
    O.(k).R = R;
    % ---- secondary: by start phase (DESCRIPTIVE; censoring caveat)
    for ph = ["light", "dark"]
        if ph == "light", m = zt < 12; else, m = zt >= 12; end
        if nnz(m) < 50, fprintf('%s/%s: %d bouts, skipped\n', k, ph, nnz(m)); continue; end
        Hp = BoutHazard(x(m), xmin, SamplingInterval=1, NumBins=16, MinAtRisk=10, NullSurvivalHandle=nullS, ...
                        NullReplicates=999, RandomSeed=seed);
        fprintf('[secondary, descriptive] %s bouts starting in %s: n %d, median %g s, RiseNullP(16 bins, pooled null) %.4f\n', ...
                k, ph, nnz(m), median(x(m)), Hp.RiseNullP);
        O.(k).(char("H16_" + ph)) = Hp;
    end
end
save(fullfile(outDir, 'o1a_results.mat'), 'O', '-v7.3');
fprintf('\nsaved %s\n', fullfile(outDir, 'o1a_results.mat'));
end
