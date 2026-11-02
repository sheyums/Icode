function X = sqwBatchSymbols(input, options)
%SQWBATCHSYMBOLS  1-minute S/Q/W symbols for every fly in a folder, written in the lab's coding (S = 0, Q = 1, W = 2).
%
%   X = sqwBatchSymbols(folder)
%   X = sqwBatchSymbols(folder, LightDark=false, Theta=0.5, OutDir=..., Name="symbols")
%
%   input      a folder of 1-s count files (*.txt, searched recursively, SQWLISTFILES) or a
%              list of such files; each file one count per line starting at ZT0 (READCOUNTS)
%   LightDark  true (LD record): Q/W reference r_f per light and per dark block, phases
%              1 = light (ZT 0-12), 2 = dark; false (DD record): one r_f, one phase (1)
%   Theta      a wake minute is Q if its active fraction < Theta * r_f (default 0.5)
%   OutDir     where to write; default "<folder>_symbols" BESIDE the input folder (not inside
%              it, so a later run on the folder does not read the outputs back)
%
%   Scoring is SQWSYMBOLS exactly as every ledger number used it: a second is active if its
%   count > 0; sleep = inactive runs >= 5 min with the lab's single-second bridging; a minute
%   is S if >= 30 of its seconds are sleep; whole days only; a wake minute in a phase with no
%   reference (fewer than 3 wake events) is missing. NOTE this is not SLEEPTRACE + IRONOUT
%   (which bridges any wake run up to maxWakeSeconds), so bout counts can differ from those.
%
%   Writes, per fly, <file>_symbols.csv with columns minute (1-based), zt_h (clock time of the
%   minute's start, h), phase, symbol (0 = S, 1 = Q, 2 = W, NaN = missing), active_frac
%   (fraction of the minute's seconds with count > 0) and sleep_s (seconds of the minute in a
%   sleep bout); and <Name>_summary.csv, one row per fly: minutes, days, counts of S/Q/W/missing
%   overall and per phase, r_f per phase, number of sleep bouts, Theta.
%   X  struct: summary (table), symbols (cell, lab coding), phase (cell), files, outDir.
%
%   See also SQWBATCHHMM, SQWSYMBOLS, SQWTOLABCODE, LABCODETOSQW.

arguments
    input
    options.LightDark (1,1) logical = true
    options.Theta (1,1) double {mustBePositive} = 0.5
    options.OutDir (1,1) string = ""
    options.Name (1,1) string = "symbols"
end
files = sqwListFiles(input);
outDir = options.OutDir;
if strlength(outDir) == 0
    if ~(isscalar(string(input)) && isfolder(input)), error('sqw:outDir', 'Give OutDir when the input is a file list.'); end
    in = char(input); if any(in(end) == '/\'), in = in(1:end - 1); end
    outDir = string(in) + "_symbols";
end
if ~isfolder(outDir), mkdir(outDir); end
if options.LightDark, names = ["light" "dark"]; else, names = "whole"; end
nF = numel(files); rows = cell(nF, 1); X.symbols = cell(nF, 1); X.phase = cell(nF, 1);
for i = 1:nF
    [y, refs, info] = sqwSymbols(readCounts(files(i)), Theta=options.Theta, LightDark=options.LightDark);
    T = numel(y);
    if options.LightDark, ph = minutePhases(T); else, ph = minutePhases(T, Scheme="none"); end
    s = sqwToLabCode(y); m = (1:T)';
    [~, nm, ex] = fileparts(files(i)); base = erase(string(nm) + ex, ".txt");
    writetable(table(m, mod((m - 1) / 60, 24), ph, s, info.activeFraction(1:T), info.sleepSeconds(1:T), ...
        'VariableNames', {'minute', 'zt_h', 'phase', 'symbol', 'active_frac', 'sleep_s'}), fullfile(outDir, base + "_symbols.csv"));
    r = table(string(nm) + ex, T, T / 1440, nnz(s == 0), nnz(s == 1), nnz(s == 2), nnz(isnan(s)), size(info.bouts, 1), options.Theta, ...
        'VariableNames', {'file', 'minutes', 'days', 'S', 'Q', 'W', 'missing', 'sleepBouts', 'theta'});
    for p = 1:numel(names)
        q = s(ph == p);
        r.("S_" + names(p)) = nnz(q == 0); r.("Q_" + names(p)) = nnz(q == 1);
        r.("W_" + names(p)) = nnz(q == 2); r.("missing_" + names(p)) = nnz(isnan(q));
        r.("r_" + names(p)) = refs.(names(p));
    end
    rows{i} = r; X.symbols{i} = s; X.phase{i} = ph;
end
X.summary = vertcat(rows{:}); X.files = files; X.outDir = outDir;
writetable(X.summary, fullfile(outDir, options.Name + "_summary.csv"));
fprintf('sqwBatchSymbols: %d flies -> %s\n', nF, outDir);
end
