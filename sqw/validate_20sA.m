function V = validate_20sA(dataRoot, outDir)
%VALIDATE_20SA  Stage 1 on raw data: reproduce the 20s-A ledger medians (hmm_audit/RESULTS_LEDGER.md).
%
%   V = validate_20sA(dataRoot, outDir)
%
%   dataRoot  folder with one subfolder per genotype holding the count files
%             cClaude used: CS, per0LD, per0DD, fumin, sss54 (the 54 sss files
%             of its run), CSDD
%   outDir    where the CSVs and figures go
%
%   Medians are deterministic and must agree to the ledger's printed
%   precision; bootstrap intervals are not compared (numpy and MATLAB
%   streams differ). Prints every comparison and returns them in V.cmp.

arguments
    dataRoot (1,1) string {mustBeFolder}
    outDir (1,1) string
end
here = fileparts(mfilename('fullpath'));
assert(isequal(fileparts(which('sqwObservableGroup')), here), 'Run from the sqw folder.');
spec = {"CS", "CS", true; "per0LD", "per0LD", true; "per0DD", "per0DD", false; ...
        "fumin", "fumin", true; "sss", "sss54", true; "CSDD", "CSDD", false};
for k = 1:size(spec, 1)
    t0 = tic;
    G(k) = sqwObservableGroup(fullfile(dataRoot, spec{k, 2}), spec{k, 1}, LightDark=spec{k, 3}, ...
        SixHour=spec{k, 1} == "CS", OutDir=outDir); %#ok<AGROW>
    fprintf('%-6s %2d flies, %.0f s\n', spec{k, 1}, numel(G(k).files), toc(t0));
end
Oc = sqwOnsagerGroup(G, Points="count", WaldSet=["per0LD" "per0DD" "fumin" "sss"], OutDir=outDir);
Oe = sqwOnsagerGroup(G, Points="chain", WaldSet=["per0LD" "per0DD" "fumin" "sss"], OutDir=outDir);

% ledger, 20s-A A1 (median columns) and A3, in BL order
L = [ ...  n nCyc  c_obs   J      EPR    Jch    Ach   EPRch  E21    E31
    25 24  0.042  0.062  0.038  0.060  0.20  0.037 -0.491 -0.909   % CS light
    25 25  0.147  0.167  0.116  0.167  0.69  0.115 -0.582 -1.782   % CS dark
    27 26  0.042  0.062  0.012  0.049  0.26  0.018 -0.572 -1.714   % per0LD light
    27 27  0.096  0.167  0.073  0.141  0.46  0.072 -0.507 -1.584   % per0LD dark
    25 24  0.118  0.167  0.110  0.171  0.60  0.111 -0.337 -1.765   % per0DD
    49 47  0.167  0.125  0.132  0.124  1.01  0.163 -1.031  0.641   % fumin light
    53 49  0.250  0.083  0.093  0.089  1.36  0.132 -1.459  1.002   % fumin dark
    50 37 -0.012  0.000  0.016 -0.042 -0.41  0.046 -1.333  1.015   % sss light
    53 48  0.000  0.000  0.029 -0.018 -0.07  0.043 -1.127  1.367   % sss dark
    25 25  0.103  0.115  0.053  0.111  0.50  0.053 -0.730 -1.263]; % CSDD
dec = [0 0 3 3 3 3 2 3 3 3];
names = ["n" "nCycle" "c_obs" "J" "EPR" "J_chain" "A_chain" "EPR_chain" "E21" "E31"];
cols = ["n" "nCycle" "c_obs_med" "J_med" "EPR_med" "J_chain_med" "A_chain_med" "EPR_chain_med" "E21_med" "E31_med"];
blk = [1 1; 1 2; 2 1; 2 2; 3 1; 4 1; 4 2; 5 1; 5 2; 6 1];
cmp = table(); nBad = 0;
fprintf('\nA1 / A3 medians: port vs ledger (rounded to the ledger''s decimals)\n');
for b = 1:size(blk, 1)
    S = G(blk(b, 1)).summary; s = S(S.phase == blk(b, 2), :);
    v = arrayfun(@(c) s.(c), cols);
    ok = abs(v - L(b, :)) <= 10 .^ -dec / 2 + 1e-9;      % agrees at the printed precision
    nBad = nBad + nnz(~ok);
    lab = G(blk(b, 1)).genotype + " " + s.phaseName;
    mark = repmat("", 1, 10); mark(~ok) = "*";
    fprintf('%-13s port   %s\n%-13s ledger %s %s\n', lab, strjoin(arrayfun(@(a, k) string(sprintf('%.*f', k, a)), v, dec) + mark, '  '), ...
        '', strjoin(compose('%g', L(b, :)), '  '), okStr(all(ok)));
    r = array2table([v; L(b, :)], 'VariableNames', cellstr(names)); r.block = [lab; lab]; r.source = ["port"; "ledger"];
    cmp = [cmp; r]; %#ok<AGROW>
end

% C3: CS 6-h blocks
L6 = [0.018 -0.042 0.045; 0.043 0.125 0.082; 0.152 0.167 0.165; 0.133 0.167 0.096];
S6 = G(1).sixSummary;
fprintf('\nC3 CS 6-h blocks (c_obs, J, EPR medians): port vs ledger\n');
for b = 1:4
    v = [S6.c_obs_med(b), S6.J_med(b), S6.EPR_med(b)];
    ok = abs(round(v, 3) - L6(b, :)) <= 5e-4 + 1e-12; nBad = nBad + nnz(~ok);
    fprintf('  %-8s port %s  ledger %s %s\n', S6.phaseName(b), mat2str(round(v, 3)), mat2str(L6(b, :)), okStr(all(ok)));
end

% A2: Onsager L (count and chain points), n
LA2 = [49 0.272 49 0.281; 53 0.286 53 0.289; 24 0.352 24 0.352; 96 0.243 84 0.232; 85 0.145 64 0.145; 25 0.262 25 0.265];
fprintf('\nA2 Onsager L: port (count n, L | chain n, L) vs ledger\n');
for g = 1:6
    v = [Oc.n(g), Oc.L(g), Oe.n(g), Oe.L(g)];
    ok = v([1 3]) == LA2(g, [1 3]) & abs(round(v([2 4]), 3) - LA2(g, [2 4])) <= 5e-4 + 1e-12;
    nBad = nBad + nnz(~ok);
    fprintf('  %-6s port %3d %.3f | %3d %.3f   ledger %3d %.3f | %3d %.3f %s\n', Oc.genotype(g), v, LA2(g, :), okStr(all(ok)));
end
fprintf('\n%d disagreement(s) with the ledger medians\n', nBad);
V = struct('G', G, 'onsagerCount', Oc, 'onsagerChain', Oe, 'cmp', cmp, 'nDisagree', nBad);
save(fullfile(outDir, 'validate_20sA.mat'), 'V', '-v7.3');
try
    f1 = plotSqwGroup(G); exportgraphics(f1, fullfile(outDir, 'observable_by_genotype.png'), 'Resolution', 150);
    f2 = plotSqwOnsager(Oc); exportgraphics(f2, fullfile(outDir, 'onsager_count.png'), 'Resolution', 150);
catch err
    fprintf('figures not written: %s\n', err.message);
end
end

function s = okStr(ok)
if ok, s = 'ok'; else, s = '<-- DIFFERS'; end
end
