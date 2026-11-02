function [y, ph, nP, refs, names] = sqwFlyInput(file, lightDark)
%SQWFLYINPUT  Symbols and phases of one fly exactly as the 20b/20p/20q analyses read them.
%
%   [y, ph, nP, refs, names] = sqwFlyInput(file, lightDark)
%
%   lightDark true: LD record, Q/W reference per light/dark block, phases 12-h
%             light (1) / dark (2); false: DD record, one whole-record
%             reference, one phase
%   y      symbols (SQWSYMBOLS, theta 0.5, whole days) in the port's INTERNAL coding 1 = Q, 2 = W,
%          3 = S, 0 = missing; SQWTOLABCODE(y) gives the lab coding S 0, Q 1, W 2, NaN
%   refs   struct of Q/W references (NaN: none, so the phase is excluded)
%   names  (nP x 1) phase names: "light"; "dark" or "whole"
%
%   See also SQWSYMBOLS, MINUTEPHASES.

[y, refs] = sqwSymbols(readCounts(file), LightDark=lightDark);
if lightDark
    ph = minutePhases(numel(y)); nP = 2; names = ["light"; "dark"];
else
    ph = minutePhases(numel(y), Scheme="none"); nP = 1; names = "whole";
end
end
