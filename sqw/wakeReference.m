function r = wakeReference(ev, options)
%WAKEREFERENCE  Fly reference activity r_f for one block (sqw_ref.wake_reference).
%
%   r = wakeReference(ev)               % whole record
%   r = wakeReference(ev, Phase=p)      % p = 1 light, 2 dark
%   r = wakeReference(ev, Phase=p, MinEvents=3)
%
%   r is the median of ev.ref over the events whose waking OR onset falls in
%   the phase; NaN if fewer than MinEvents events qualify. A fly-phase with
%   NaN has no Q/W reference: its wake minutes cannot be called Q or W.
%
%   See also WAKEBOUTEVENTS, SQWSYMBOLS.

arguments
    ev (:,1) struct
    options.Phase double {mustBeScalarOrEmpty} = []
    options.MinEvents (1,1) double {mustBePositive, mustBeInteger} = 3
end

if isempty(ev)
    r = NaN; return
end
ref = [ev.ref]';
if isempty(options.Phase)
    use = true(size(ref));
else
    use = [ev.wakePhase]' == options.Phase | [ev.onsetPhase]' == options.Phase;
end
if nnz(use) >= options.MinEvents
    r = median(ref(use));
else
    r = NaN;
end
end
