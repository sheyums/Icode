function bouts = sleepBouts(x, options)
%SLEEPBOUTS  Sleep bouts from 1-s counts, with the lab's single-second bridging (sqw_ref.sleep_bouts).
%
%   bouts = sleepBouts(x)
%   bouts = sleepBouts(x, Bridge=true, MinLen=300)
%
%   x       1-s movement counts (a second is ACTIVE if its count > 0)
%   Bridge  reset an isolated active second to inactive when the inactive runs
%           on BOTH sides last >= MinLen s (default true)
%   MinLen  shortest sleep bout, seconds (default 300, the 5-minute criterion)
%
%   bouts   (nBouts x 2) [start, length], start 1-BASED (second 1 = the file's
%           first line); the golden files store 0-based starts
%
%   Bridging is lab convention (b): flank lengths are measured on the
%   UNBRIDGED trace, so 300+1+300+1+300 becomes one 902-s bout, while
%   301+1+100 and 300+2+300 are not bridged. A run touching the start or end
%   of the trace counts as a flank. A bout touching either edge of the trace
%   is KEPT (it is censored, but sqw_ref keeps it and every downstream count
%   was made that way).
%
%   See also WAKEBOUTEVENTS, SQWSYMBOLS.

arguments
    x (:,1) double {mustBeNonnegative}
    options.Bridge (1,1) logical = true
    options.MinLen (1,1) double {mustBePositive, mustBeInteger} = 300
end

active = x > 0;
if options.Bridge && numel(x) >= 3
    [st, len, isRest] = runsOf(~active);
    i = (2:numel(st) - 1)';                         % candidate centres
    hit = i(~isRest(i) & len(i) == 1 & isRest(i - 1) & isRest(i + 1) & ...
        len(i - 1) >= options.MinLen & len(i + 1) >= options.MinLen);
    active(st(hit)) = false;
end
[st, len, isRest] = runsOf(~active);
keep = isRest & len >= options.MinLen;
bouts = [st(keep), len(keep)];
if isempty(bouts), bouts = zeros(0, 2); end
end

function [st, len, val] = runsOf(v)
% Maximal runs of a logical vector: 1-based starts, lengths, values.
if isempty(v)
    st = zeros(0, 1); len = st; val = false(0, 1); return
end
st = [1; find(diff(v) ~= 0) + 1];
len = diff([st; numel(v) + 1]);
val = v(st);
end
