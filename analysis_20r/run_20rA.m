function A = run_20rA(stageDir)
% RUN_20RA  20r-A on CS DD (PRESPEC 2026-10-31 20r-A, with the extractor of
% "20r amendment: stored per0 DD bouts lost"). Run BY NAME from the Icode
% folder, with analysis_20r on the path, AFTER the per0 DD replication.
%
%   stageDir : folder holding ONLY the 25 CS DD files, cut to 345600 samples.
%
% Steps, as pre-specified: (1) rise test on wake bouts at 12/16/20/24 bins;
% (2) per-fly influence, 2026-09-15 design; (3) the same rise test on SLEEP
% bouts (control); (4) CompareBoutModels(x, 2, SamplingInterval=1,
% MaxComponents=5), default library. Output: analysis_20r/CSDD_bouts_20r.mat,
% analysis_20r/run_20rA_results.mat, and the console log via diary.

assert(isequal(fileparts(which('BoutHazard')), pwd), ...
    'Run from the Icode folder: BoutHazard resolves to %s', which('BoutHazard'));
assert(isequal(fileparts(which('CompareBoutModels')), pwd), ...
    'CompareBoutModels resolves to %s', which('CompareBoutModels'));
here = fileparts(mfilename('fullpath'));

B = extract_20r(stageDir, fullfile(here, 'CSDD_bouts_20r.mat'));
nW = nnz(B.Wake >= B.xminWake); nS = nnz(B.Sleep >= B.xminSleep);
fprintf('CS DD: %d flies, %d wake bouts, %d sleep bouts, %.2f fly-days of record\n', ...
    numel(B.files), nW, nS, sum(B.traceLen) / 86400);

A = struct();
A.wake      = rise_test_20r(B.Wake,  B.xminWake,  'CS DD wake (20r-A)');
A.influence = influence_20r(B.Wake,  B.xminWake,  'CS DD wake (20r-A)');
A.sleep     = rise_test_20r(B.Sleep, B.xminSleep, 'CS DD sleep, control (20r-A)');

x = B.Wake(isfinite(B.Wake) & B.Wake >= B.xminWake);
R = CompareBoutModels(x, B.xminWake, SamplingInterval=1, MaxComponents=5, Plot=false);
fprintf('\nCompareBoutModels, CS DD wake (read in CLAUDE.md order)\n');
fprintf('R.Skipped:\n'); disp(R.Skipped);
fprintf('R.Table:\n'); disp(R.Table);
fprintf('R.Selected:\n'); disp(R.Selected);
fprintf('R.SelectionPath:\n'); disp(R.SelectionPath);
fprintf('R.GoF:\n'); disp(R.GoF);
if isfield(R, 'SelectedFit') && ~isempty(R.SelectedFit)
    fprintf('R.SelectedFit:\n'); disp(R.SelectedFit);
    if isfield(R.SelectedFit, 'Diagnostics'), disp(R.SelectedFit.Diagnostics); end
end
A.library = rmfield(R, intersect(fieldnames(R), {'Figure'}));
save(fullfile(here, 'run_20rA_results.mat'), 'A', '-v7');
end
