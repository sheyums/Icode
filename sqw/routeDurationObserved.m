function R = routeDurationObserved(x, options)
%ROUTEDURATIONOBSERVED  Entry and exit route of every observed sleep bout, with its length (sqw_ref.route_duration_observed; PRESPEC 20k).
%
%   R = routeDurationObserved(x)
%   R = routeDurationObserved(x, Theta=0.5, LightDark=true, MinWake=600)
%
%   x  1-s counts from ZT0 (READCOUNTS), the whole record
%
%   R  (nBouts x 1) struct, one per bridged sleep bout (SLEEPBOUTS) with a sleep bout before
%      AND after it, skipped if its onset phase has no Q/W reference:
%        phase          1 light / 2 dark from the onset second's clock time (1 for DD)
%        sleep_min, wake_before_min, wake_after_min   lengths in minutes
%        q_in, q_out    1 = quiet route, 0 = active, NaN = not classified: the wake bout
%                       before (after) is shorter than MinWake s. Quiet if the active fraction
%                       in seconds 61-300 before onset (after waking) is < Theta r_f, with r_f
%                       the reference of the ONSET's phase (WAKEREFERENCE).
%   Windows in 0-based seconds as sqw_ref: before [s - 300, s - 60), after [e + 60, e + 300),
%   s = onset, e = first second after the bout.
%
%   See also ROUTEDURATIOND, SLEEPBOUTS, WAKEBOUTEVENTS, WAKEREFERENCE.

arguments
    x (:,1) double
    options.Theta (1,1) double {mustBePositive} = 0.5
    options.LightDark (1,1) logical = true
    options.MinWake (1,1) double {mustBePositive} = 600
end
bouts = sleepBouts(x); mv = double(x > 0); ev = wakeBoutEvents(x, bouts);
if options.LightDark, ref = [wakeReference(ev, Phase=1), wakeReference(ev, Phase=2)];
else, ref = wakeReference(ev); end
nb = size(bouts, 1); c = cell(max(nb - 2, 0), 1);
for k = 2:nb - 1
    S = bouts(k, 1); n = bouts(k, 2); s0 = S - 1; e0 = s0 + n;            % 0-based onset and end
    wb = s0 - (bouts(k - 1, 1) - 1 + bouts(k - 1, 2)); wa = (bouts(k + 1, 1) - 1) - e0;
    if options.LightDark, p = 1 + (mod(s0, 86400) >= 43200); else, p = 1; end
    r = ref(p);
    if ~isfinite(r), continue, end
    qi = NaN; qo = NaN;
    if wb >= options.MinWake, qi = double(mean(mv(s0 - 300 + 1:s0 - 60)) < options.Theta * r); end
    if wa >= options.MinWake, qo = double(mean(mv(e0 + 60 + 1:e0 + 300)) < options.Theta * r); end
    c{k - 1} = struct('phase', p, 'sleep_min', n / 60, 'wake_before_min', wb / 60, 'wake_after_min', wa / 60, ...
        'q_in', qi, 'q_out', qo);
end
c = c(~cellfun(@isempty, c));
if isempty(c)
    R = struct('phase', {}, 'sleep_min', {}, 'wake_before_min', {}, 'wake_after_min', {}, 'q_in', {}, 'q_out', {});
else
    R = vertcat(c{:});
end
end
