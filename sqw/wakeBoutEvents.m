function ev = wakeBoutEvents(x, bouts, options)
%WAKEBOUTEVENTS  Activity of long wake bouts between two sleep bouts (sqw_ref.wake_bout_events; PRESPEC 19).
%
%   ev = wakeBoutEvents(x, bouts)
%   ev = wakeBoutEvents(x, bouts, MinWake=600, RefHalfWidth=150, ...
%                       TransitWindow=[60 300], DaySeconds=86400, LightSeconds=43200)
%
%   x       1-s counts (a second is active if its count > 0)
%   bouts   SLEEPBOUTS output, [start (1-based), length]
%
%   ev      (nEvents x 1) struct, one per wake bout of >= MinWake s lying
%           between two consecutive sleep bouts:
%             ref         active fraction in the bout's middle 2*RefHalfWidth s
%             wakePhase   1 light / 2 dark, clock phase of the bout's first second
%             wakeAct     active fraction in seconds 61-300 after waking
%             onsetPhase  1 / 2, clock phase of the next sleep onset
%             onsetAct    active fraction in seconds 61-300 before that onset
%             wakeStart, wakeEnd   first and last second of the wake bout (1-based)
%
%   Windows, in 0-based seconds as in sqw_ref, with w0 = first wake second and
%   w1 = first second of the next sleep bout:
%     ref    [mid - 150, mid + 150),   mid = floor((w0 + w1)/2)
%     exit   [w0 + 60, w0 + 300)       approach  [w1 - 300, w1 - 60)
%   A second s (0-based) is DARK if mod(s, DaySeconds) >= LightSeconds.
%
%   See also WAKEREFERENCE, SQWSYMBOLS.

arguments
    x (:,1) double
    bouts (:,2) double
    options.MinWake (1,1) double {mustBePositive} = 600
    options.RefHalfWidth (1,1) double {mustBePositive, mustBeInteger} = 150
    options.TransitWindow (1,2) double {mustBeNonnegative, mustBeInteger} = [60 300]
    options.DaySeconds (1,1) double {mustBePositive} = 86400
    options.LightSeconds (1,1) double {mustBePositive} = 43200
end

tw = options.TransitWindow;
if tw(2) <= tw(1) || 2 * options.RefHalfWidth > options.MinWake || tw(2) > options.MinWake
    error('wakeBoutEvents:windows', ...
        'Windows must fit inside the shortest wake bout: TransitWindow(2) and 2*RefHalfWidth <= MinWake.');
end
mv = double(x > 0);
h = options.RefHalfWidth;
dark = @(s0) 1 + double(mod(s0, options.DaySeconds) >= options.LightSeconds);
frac = @(lo0, hi0) mean(mv(lo0 + 1:hi0));          % mean over 0-based [lo0, hi0)

ev = repmat(struct('ref', NaN, 'wakePhase', NaN, 'wakeAct', NaN, 'onsetPhase', NaN, ...
    'onsetAct', NaN, 'wakeStart', NaN, 'wakeEnd', NaN), 0, 1);
for k = 1:size(bouts, 1) - 1
    w0 = bouts(k, 1) - 1 + bouts(k, 2);            % 0-based first wake second
    w1 = bouts(k + 1, 1) - 1;                      % 0-based next onset
    if w1 - w0 < options.MinWake, continue, end
    mid = floor((w0 + w1) / 2);
    e.ref = frac(mid - h, mid + h);
    e.wakePhase = dark(w0);
    e.wakeAct = frac(w0 + tw(1), w0 + tw(2));
    e.onsetPhase = dark(w1);
    e.onsetAct = frac(w1 - tw(2), w1 - tw(1));
    e.wakeStart = w0 + 1;
    e.wakeEnd = w1;
    ev(end + 1, 1) = e; %#ok<AGROW>
end
end
