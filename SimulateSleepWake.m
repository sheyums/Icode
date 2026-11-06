function fly = SimulateSleepWake(mode, seed, opts)
%SIMULATESLEEPWAKE  Simulated alternating sleep/wake records with known structure.
%
%   fly = SimulateSleepWake(mode, seed)
%   fly = SimulateSleepWake(mode, seed, Days=5, Gain=2, ...)
%
% A test and calibration utility for the homeostat analyses (BoutCrossType,
% SleepAmountACF, ProcessSProxy). It is NOT a model of flies: each mode
% builds in exactly one kind of structure, so a test's behaviour on it is
% known in advance.
%
% Every bout is exponential above its floor (sleep >= SleepMin, wake >=
% WakeMin, integer seconds), with a mean that depends on the mode:
%   "iid"        constant means: no structure at all.
%   "shared"     one latent AR(1) z (per sleep-wake cycle, coefficient Phi)
%                multiplies BOTH means by exp(Sigma*z): a slowly varying
%                state that lengthens sleep and wake together. No homeostat.
%   "homeostat"  a Borbely-type pressure S in [0, 1], rising during wake
%                (time constant Tw) and falling during sleep (Ts); the sleep
%                mean is multiplied by exp(Gain*(S - 0.5)) and the wake mean
%                by exp(-Gain*(S - 0.5)), with S taken at the bout's onset.
%   "ld"         a 12:12 time-of-day profile, no coupling: sleep mean 15 min
%                by day (ZT 0-12), rising 20 -> 80 min through the night;
%                wake mean 25 min by day, 8 min by night (as in the
%                BoutIrreversibility LD-like tests).
%
% OUTPUT  struct with column vectors
%   Dur      bout durations, seconds (integers)
%   IsSleep  true for sleep bouts
%   Onset    seconds from the record start; contiguous: Onset(i+1) =
%            Onset(i) + Dur(i)
%   Strata   floor(ZT of onset), ZT = mod(Onset/3600, 24) (record starts at
%            ZT0)
%   Truth    the mode and parameters used
%
% Uses a private RandStream ('mt19937ar', Seed); the global stream is not
% touched.
%
% CHANGE LOG
%   2026-11-06  First version.

arguments
    mode (1,1) string {mustBeMember(mode, ["iid","shared","homeostat","ld"])}
    seed (1,1) double {mustBeInteger, mustBeNonnegative}
    opts.Days (1,1) double {mustBePositive} = 5
    opts.SleepMeanMin (1,1) double {mustBePositive} = 15
    opts.WakeMeanMin (1,1) double {mustBePositive} = 15
    opts.SleepMin (1,1) double = 300
    opts.WakeMin (1,1) double = 2
    opts.Phi (1,1) double = 0.9
    opts.Sigma (1,1) double = 0.5
    opts.Tw (1,1) double {mustBePositive} = 2       % hours
    opts.Ts (1,1) double {mustBePositive} = 1       % hours
    opts.Gain (1,1) double = 2
end
rs = RandStream('mt19937ar', 'Seed', seed);
T = opts.Days * 86400;
nMax = ceil(T / 60) + 10;
dur = zeros(nMax, 1); isS = false(nMax, 1); on = zeros(nMax, 1);
t = 0; k = 0;
sleepNow = rs.rand < 0.5;
z = rs.randn; S = 0.5;
while t < T
    k = k + 1;
    zt = mod(t / 3600, 24);
    switch mode
        case "iid"
            ms = opts.SleepMeanMin; mw = opts.WakeMeanMin;
        case "shared"
            if sleepNow     % one z per cycle, updated at each sleep onset
                z = opts.Phi * z + sqrt(1 - opts.Phi^2) * rs.randn;
            end
            ms = opts.SleepMeanMin * exp(opts.Sigma * z);
            mw = opts.WakeMeanMin * exp(opts.Sigma * z);
        case "homeostat"
            ms = opts.SleepMeanMin * exp(opts.Gain * (S - 0.5));
            mw = opts.WakeMeanMin * exp(-opts.Gain * (S - 0.5));
        case "ld"
            if zt < 12
                ms = 15; mw = 25;
            else
                ms = 20 + 60 * (zt - 12) / 12; mw = 8;
            end
    end
    if sleepNow
        d = opts.SleepMin + round(60 * ms * -log(rs.rand));
    else
        d = opts.WakeMin + round(60 * mw * -log(rs.rand));
    end
    dur(k) = d; isS(k) = sleepNow; on(k) = t;
    if mode == "homeostat"
        if sleepNow
            S = S * exp(-d / (3600 * opts.Ts));
        else
            S = 1 - (1 - S) * exp(-d / (3600 * opts.Tw));
        end
    end
    t = t + d;
    sleepNow = ~sleepNow;
end
fly.Dur = dur(1:k); fly.IsSleep = isS(1:k); fly.Onset = on(1:k);
fly.Strata = floor(mod(fly.Onset / 3600, 24));
fly.Truth = opts; fly.Truth.Mode = mode;
end
