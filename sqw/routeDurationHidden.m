function R = routeDurationHidden(z, ph, options)
%ROUTEDURATIONHIDDEN  Route-vs-duration records on one hidden path (sqw_ref.route_duration_hidden; PRESPEC 20k).
%
%   R = routeDurationHidden(z, ph)
%   R = routeDurationHidden(z, ph, MinSleep=5, MinGap=10)
%
%   z    (T x 1) hidden path, 1 = H1 (quiet), 2 = H2 (active), 3 = H3 (sleep-like)
%   ph   (T x 1) phase of each minute
%   R    struct array as ROUTEDURATIONOBSERVED, one per H3 run of >= MinSleep minutes with an H3
%        run before and after it: phase = phase of its first minute, sleep_min = its length,
%        wake_before_min / wake_after_min = the gaps to the neighbouring H3 runs, and q_in / q_out
%        = 1 if the run before / after it is H1, 0 if not, NaN if that gap is < MinGap minutes.
%
%   See also ROUTEDURATIONOBSERVED, ROUTEDURATIOND.

arguments
    z (:,1) double
    ph (:,1) double
    options.MinSleep (1,1) double = 5
    options.MinGap (1,1) double = 10
end
ch = [1; find(diff(z) ~= 0) + 1; numel(z) + 1]; st = z(ch(1:end - 1)); ln = diff(ch);
si = find(st == 3); c = cell(max(numel(si) - 2, 0), 1);
for a = 2:numel(si) - 1
    k = si(a);
    if ln(k) < options.MinSleep, continue, end
    gb = ch(k) - ch(si(a - 1) + 1); ga = ch(si(a + 1)) - ch(k + 1);
    qi = NaN; qo = NaN;
    if gb >= options.MinGap, qi = double(st(k - 1) == 1); end
    if ga >= options.MinGap, qo = double(st(k + 1) == 1); end
    c{a - 1} = struct('phase', ph(ch(k)), 'sleep_min', ln(k), 'wake_before_min', gb, 'wake_after_min', ga, 'q_in', qi, 'q_out', qo);
end
c = c(~cellfun(@isempty, c));
if isempty(c)
    R = struct('phase', {}, 'sleep_min', {}, 'wake_before_min', {}, 'wake_after_min', {}, 'q_in', {}, 'q_out', {});
else
    R = vertcat(c{:});
end
end
