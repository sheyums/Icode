function S = paperFigureStyle()
%PAPERFIGURESTYLE  Common style of the Plan A figures (FIGURE_RECIPES section K).
%
%   S = paperFigureStyle()
%
%   S.order   genotype order rule: CS -> CS DD -> per0 LD -> per0 DD -> fumin -> sss
%   S.dd      true for the DD genotypes (one phase, "whole")
%   S.col     (6 x 3) RGB colours: CS #2a78d6, CS DD #6250d6, per0 LD #1baf7a,
%             per0 DD #eda100, fumin #eb6834, sss #e87ba4
%   S.name    display names (TeX interpreter)
%   S.blocks  table of the strip blocks: genotype index, phase (1 light / dark /
%             whole, 2 dark), filled (dark or DD), label; light before dark
%   S.font, S.size   font (STIX if installed, else Times New Roman), 8 pt
%   S.NBoot, S.Level  fly-bootstrap CI of the median: 2000 resamples, 95%
%
%   See also NEWPAPERFIGURE, EXPORTDIRECTIONFLUXFIGURE.

S.order = ["CS" "CSDD" "per0LD" "per0DD" "fumin" "sss"];
S.dd = [false true false true false false];
hex = ["2a78d6" "6250d6" "1baf7a" "eda100" "eb6834" "e87ba4"];
S.col = zeros(6, 3);
for k = 1:6, S.col(k, :) = sscanf(hex(k), '%2x%2x%2x')' / 255; end
S.name = ["CS" "CS DD" "\itper^0\rm LD" "\itper^0\rm DD" "\itfmn\rm" "\itsss\rm"];
g = []; p = []; f = []; lab = strings(0, 1);
for k = 1:6
    if S.dd(k)
        g(end + 1, 1) = k; p(end + 1, 1) = 1; f(end + 1, 1) = true; lab(end + 1, 1) = S.name(k); %#ok<AGROW>
    else
        g(end + 1, 1) = k; p(end + 1, 1) = 1; f(end + 1, 1) = false; lab(end + 1, 1) = S.name(k) + " light"; %#ok<AGROW>
        g(end + 1, 1) = k; p(end + 1, 1) = 2; f(end + 1, 1) = true; lab(end + 1, 1) = S.name(k) + " dark"; %#ok<AGROW>
    end
end
S.blocks = table(g, p, logical(f), lab, 'VariableNames', {'g', 'phase', 'filled', 'label'});
fonts = string(listfonts);
cand = ["STIX Two Text" "STIXGeneral" "STIX" "Times New Roman"];
S.font = cand(find(ismember(cand, fonts), 1));
if isempty(S.font), S.font = "Helvetica"; end
S.size = 8; S.NBoot = 2000; S.Level = 0.95;
end
