function R = turEstimate(y, ph, nP, options)
%TURESTIMATE  Counted cycle current and its TUR bound from the data, per phase (sqw_ref.tur_estimate; PRESPEC 20q).
%
%   R = turEstimate(y, ph, nP)
%   R = turEstimate(y, ph, nP, Taus=[180 120 360])
%
%   y     (T x 1) symbols 1 = Q, 2 = W, 3 = S, 0 = missing, 4 = unknown wake
%   ph    (T x 1) phase of each minute, 1..nP (arrival-minute convention)
%   Taus  block lengths in minutes (180 = 3 h is the 20q primary)
%
%   R  (nP x 1) struct array:
%        J     counted forward-cycle increments per hour: (sum of increments
%              in phase p) / (valid minutes in phase p / 60)
%        Taus  (1 x nTau) block lengths
%        D     (1 x nTau) Var(N_tau) / (2 tau), per hour; NaN with < 4 blocks
%        B     (1 x nTau) TUR bound J^2 / D (k_B per hour; 20q amendment C);
%              NaN unless D is finite and > 0
%        n     (1 x nTau) number of blocks used
%
%   Increments (CYCLEINCREMENTS) sit on transitions between two consecutive
%   observed minutes (not missing, not unknown wake), attributed to the
%   arrival minute. Blocks: consecutive, non-overlapping tau-minute blocks
%   tiling each run of constant phase from its start, kept only if every
%   minute's transition is valid (so a block never contains or follows a gap).
%   Var uses N - 1.
%
%   See also TUREXACT, TUREXACTCT, CYCLEINCREMENTS.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    options.Taus (1,:) double {mustBePositive, mustBeInteger} = [180 120 360]
end
T = numel(y);
if numel(ph) ~= T, error('sqw:length', 'y and ph must have the same length.'); end
dC = cycleIncrements();
obs = y >= 1 & y <= 3;
val = false(T, 1); val(2:end) = obs(2:end) & obs(1:end - 1);
inc = zeros(T, 1);
t = find(val);
inc(t) = dC(sub2ind([3 3], y(t - 1), y(t)));
cv = [0; cumsum(val)]; ci = [0; cumsum(inc)];       % prefix sums for block checks and sums
seg = find(diff(ph) ~= 0) + 1; starts = [1; seg]; ends = [seg - 1; T];
nT = numel(options.Taus);
R = repmat(struct('J', NaN, 'Taus', options.Taus, 'D', NaN(1, nT), 'B', NaN(1, nT), 'n', zeros(1, nT)), nP, 1);
for p = 1:nP
    m = ph == p & val; hours = nnz(m) / 60;
    if hours > 0, J = sum(inc(m)) / hours; else, J = NaN; end
    R(p).J = J;
    for q = 1:nT
        tau = options.Taus(q); Nb = [];
        for s = find(ph(starts) == p)'
            b0 = (starts(s):tau:ends(s) - tau + 1)';
            b0 = b0(b0 ~= 1);
            full = cv(b0 + tau) - cv(b0) == tau;     % every transition in b0..b0+tau-1 valid
            Nb = [Nb; ci(b0(full) + tau) - ci(b0(full))]; %#ok<AGROW>
        end
        R(p).n(q) = numel(Nb);
        if numel(Nb) >= 4
            D = var(Nb) / (2 * tau / 60);
            R(p).D(q) = D;
            if isfinite(D) && D > 0, R(p).B(q) = J^2 / D; end
        end
    end
end
end
