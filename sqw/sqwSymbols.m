function [y, refs, info] = sqwSymbols(x, options)
%SQWSYMBOLS  1-minute S/Q/W symbols from 1-s counts (sqw_ref.sqw_symbols; PRESPEC 20 section 1).
%
%   [y, refs, info] = sqwSymbols(x)
%   [y, refs, info] = sqwSymbols(x, Theta=0.5, LightDark=true, FullDays=true, ...
%                                UnknownWake="missing", SleepSecondsForS=30)
%
%   x            1-s counts starting at ZT0
%   Theta        a wake minute is Q if its active fraction < Theta * r_f (strict)
%   LightDark    true: r_f per light / dark block (LD records); false: one r_f
%                for the whole record (DD records)
%   FullDays     truncate the symbols to whole days (1440 minutes)
%   UnknownWake  code for a wake minute whose block has no reference r_f:
%                  "missing" (0; what every ledger number used, the default)
%                  "unknown" (4; wake of unknown type, the recommended partial
%                  observation of PORTING_GUIDE section 4)
%   SleepSecondsForS  a minute is S if at least this many of its 60 s lie in
%                a (bridged) sleep bout (default 30)
%
%   y      (nMinutes x 1) symbols: 1 = Q, 2 = W, 3 = S, 0 = missing, 4 = wake of
%          unknown type. (sqw_ref uses Q = 0, W = 1, S = 2, -1, -2.)
%   refs   struct: fields light, dark (LightDark) or whole
%   info   struct: bouts (on the whole-minute trace, 1-based), events,
%          sleepSeconds and activeFraction per minute (before FullDays),
%          nMinutesRaw
%
%   Minute m (1-based) covers seconds 60(m-1)+1 .. 60m and is LIGHT if
%   mod(60(m-1), 86400) < 43200. The trace is cut to whole minutes BEFORE the
%   sleep bouts are found, exactly as sqw_ref does.
%
%   See also SLEEPBOUTS, WAKEBOUTEVENTS, WAKEREFERENCE, MINUTEPHASES.

arguments
    x (:,1) double {mustBeNonnegative}
    options.Theta (1,1) double {mustBePositive} = 0.5
    options.LightDark (1,1) logical = true
    options.FullDays (1,1) logical = true
    options.UnknownWake (1,1) string {mustBeMember(options.UnknownWake, ["missing" "unknown"])} = "missing"
    options.SleepSecondsForS (1,1) double {mustBePositive} = 30
end

Q = 1; W = 2; S = 3;
nm = floor(numel(x) / 60);
x = x(1:nm * 60);
bouts = sleepBouts(x, Bridge=true);
sl = false(size(x));
for k = 1:size(bouts, 1)
    sl(bouts(k, 1):bouts(k, 1) + bouts(k, 2) - 1) = true;
end
ssec = sum(reshape(sl, 60, nm), 1)';
act = mean(reshape(double(x > 0), 60, nm), 1)';
ev = wakeBoutEvents(x, bouts);
if options.LightDark
    refs = struct('light', wakeReference(ev, Phase=1), 'dark', wakeReference(ev, Phase=2));
    isLight = mod(60 * (0:nm - 1)', 86400) < 43200;
    r = refs.dark * ones(nm, 1);
    r(isLight) = refs.light;
else
    refs = struct('whole', wakeReference(ev));
    r = refs.whole * ones(nm, 1);
end
if options.UnknownWake == "missing", unk = 0; else, unk = 4; end

y = W * ones(nm, 1);
y(ssec >= options.SleepSecondsForS) = S;
wake = ssec < options.SleepSecondsForS;
y(wake & isnan(r)) = unk;
y(wake & ~isnan(r) & act < options.Theta * r) = Q;
if options.FullDays
    y = y(1:floor(nm / 1440) * 1440);
end
info = struct('bouts', bouts, 'events', ev, 'sleepSeconds', ssec, ...
    'activeFraction', act, 'nMinutesRaw', nm);
end
