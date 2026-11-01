function C = roundTripsHidden(z, ph, nP, options)
%ROUNDTRIPSHIDDEN  Round-trip types around sleep runs on one hidden path (sqw_ref.round_trips_hidden; PRESPEC 20g).
%
%   C = roundTripsHidden(z, ph, nP)
%   C = roundTripsHidden(z, ph, nP, MinSleep=5, MinGap=10)     % matched definition (20g amendment B)
%
%   z   (T x 1) path in 1..3 (H1 plays Q, H2 plays A, H3 = sleep)
%   C   (nP x 4) counts [forward reverse direct quiet_both] = (entry, exit) of
%       (Q, A), (A, Q), (A, A), (Q, Q), by phase of the H3 run's first minute
%
%   Every H3 run with a neighbouring run on both sides (default), or, with
%   MinGap > 0, H3 runs >= MinSleep minutes with >= MinGap minutes free of H3
%   before and after (the previous and next H3 runs must exist).
%
%   See also ROUNDTRIPSTATS.

arguments
    z (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    options.MinSleep (1,1) double {mustBeNonnegative} = 0
    options.MinGap (1,1) double {mustBeNonnegative} = 0
end
ch = [1; find(diff(z) ~= 0) + 1; numel(z) + 1];        % run starts, plus end sentinel
st = z(ch(1:end - 1)); ln = diff(ch);
sIdx = find(st == 3);
C = zeros(nP, 4);
for a = 1:numel(sIdx)
    k = sIdx(a);
    if k == 1 || k == numel(st) || ln(k) < options.MinSleep, continue, end
    if options.MinGap > 0
        if a == 1 || a == numel(sIdx), continue, end
        if ch(k) - ch(sIdx(a - 1) + 1) < options.MinGap || ch(sIdx(a + 1)) - ch(k + 1) < options.MinGap
            continue
        end
    end
    qi = st(k - 1) == 1; qo = st(k + 1) == 1;
    col = [3 2; 1 4];                                   % (qi+1, qo+1): (A,A) direct, (A,Q) reverse, (Q,A) forward, (Q,Q)
    C(ph(ch(k)), col(qi + 1, qo + 1)) = C(ph(ch(k)), col(qi + 1, qo + 1)) + 1;
end
end
