function sur = RenewalSurrogate(fly, rs)
%RENEWALSURROGATE  One renewal surrogate of an alternating sleep/wake record.
%
%   sur = RenewalSurrogate(fly, rs)
%
% Within each CONTIGUOUS segment of the record (bouts i, i+1 with
% Onset(i+1) == Onset(i) + Dur(i)), the sleep durations are permuted among
% the segment's sleep positions and the wake durations among its wake
% positions; the onsets are then rebuilt end to end from the segment's
% original start. The alternation, the first state, the segment's start
% and total length, and both bout distributions are kept exactly; every
% dependence between bouts is removed. This is the null of a renewal
% process with the fly's own bout distributions (SleepAmountACF,
% ProcessSProxy).
%
% fly: struct with column vectors Dur, IsSleep, Onset (seconds). Other
% fields are copied unchanged. rs: a RandStream (the caller's, so that
% different analyses can draw IDENTICAL surrogates from the same stream).
%
% CHANGE LOG
%   2026-11-06  First version.

seg = segments(fly);
d = fly.Dur;
on = fly.Onset;
for s = 1:size(seg, 1)
    ix = (seg(s, 1):seg(s, 2)).';
    for type = [true false]
        p = ix(fly.IsSleep(ix) == type);
        d(p) = d(p(randperm(rs, numel(p))));
    end
    on(ix) = fly.Onset(ix(1)) + [0; cumsum(d(ix(1:end-1)))];
end
sur = fly;
sur.Dur = d;
sur.Onset = on;
end

function seg = segments(fly)
% [first last] bout index of each contiguous run.
brk = find(fly.Onset(2:end) ~= fly.Onset(1:end-1) + fly.Dur(1:end-1));
st = [1; brk + 1];
en = [brk; numel(fly.Dur)];
seg = [st en];
end
