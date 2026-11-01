function c = cObsMinutes(y, ph, nP)
%COBSMINUTES  Observable cycle current on the minute symbols (sqw_ref.c_obs_minutes; PRESPEC 20d).
%
%   c = cObsMinutes(y, ph, nP)
%
%   c(p) = P(last wake minute before a sleep onset is Q)
%          - P(first wake minute after a waking is Q),
%   over the onsets and wakings whose SECOND minute (the S minute of an
%   onset, the wake minute of a waking) lies in phase p. A "wake minute" is
%   Q or W; missing (0) and unknown-wake (4) minutes are neither, so a
%   transition touching one is skipped. NaN if phase p has no onset or no
%   waking. Positive = forward cycle W -> Q -> S -> W.
%
%   See also OBSERVABLEPAIRCOUNTS, OBSERVABLEQUANTITIES.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
end
if numel(y) ~= numel(ph), error('sqw:length', 'y and ph must have the same length.'); end

Q = 1; S = 3;
isWake = @(v) v == 1 | v == 2;
t = (2:numel(y))';
on = t(isWake(y(t - 1)) & y(t) == S);
wk = t(y(t - 1) == S & isWake(y(t)));
c = NaN(nP, 1);
for p = 1:nP
    o = on(ph(on) == p); w = wk(ph(wk) == p);
    if ~isempty(o) && ~isempty(w)
        c(p) = mean(y(o - 1) == Q) - mean(y(w) == Q);
    end
end
end
