function D = routeDurationD(R, relation, phase, options)
%ROUTEDURATIOND  Route-vs-duration statistic D of one fly and phase (sqw_ref.route_duration_D; PRESPEC 20k).
%
%   D = routeDurationD(R, relation, phase)
%   D = routeDurationD(R, relation, phase, MinEach=2)
%
%   R         ROUTEDURATIONOBSERVED records
%   relation  "R1_exit_vs_sleep" (q_out, sleep_min), "R2_entry_vs_wake_before" (q_in,
%             wake_before_min) or "R3_entry_vs_sleep" (q_in, sleep_min)
%   phase     1 light / 2 dark (1 for DD), the onset's phase
%
%   D = mean ln(length | quiet route) - mean ln(length | active route), over the records of
%   that phase with a classified route; NaN if fewer than MinEach of either kind.
%   D < 0: quiet routes go with shorter bouts.
%
%   See also ROUTEDURATIONOBSERVED.

arguments
    R struct
    relation (1,1) string {mustBeMember(relation, ["R1_exit_vs_sleep" "R2_entry_vs_wake_before" "R3_entry_vs_sleep"])}
    phase (1,1) double {mustBePositive, mustBeInteger}
    options.MinEach (1,1) double {mustBePositive, mustBeInteger} = 2
end
[lab, dur] = routeFields(relation);
if isempty(R), D = NaN; return, end
q = [R.(lab)]'; t = [R.(dur)]'; k = [R.phase]' == phase & ~isnan(q);
q = q(k); t = log(t(k));
if nnz(q == 1) >= options.MinEach && nnz(q == 0) >= options.MinEach
    D = mean(t(q == 1)) - mean(t(q == 0));
else
    D = NaN;
end
end

function [lab, dur] = routeFields(relation)
switch relation
    case "R1_exit_vs_sleep", lab = 'q_out'; dur = 'sleep_min';
    case "R2_entry_vs_wake_before", lab = 'q_in'; dur = 'wake_before_min';
    otherwise, lab = 'q_in'; dur = 'sleep_min';
end
end
