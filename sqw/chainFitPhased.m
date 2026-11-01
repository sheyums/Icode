function [ll, pi1, P, info] = chainFitPhased(y, ph, nP, options)
%CHAINFITPHASED  First-order Markov chain on the observed S/Q/W symbols, per phase (sqw_ref.mc_fit; PRESPEC 20n).
%
%   [ll, pi1, P, info] = chainFitPhased(y, ph, nP)
%   [ll, pi1, P, info] = chainFitPhased(y, ph, nP, Tol=1e-10, MaxIter=500)
%
%   y    (T x 1) symbols 1 = Q, 2 = W, 3 = S, 0 = missing, 4 = unknown wake
%   ph   (T x 1) phase of each minute, 1..nP (arrival-minute convention)
%
%   ll   log-likelihood of the fitted chain
%   pi1  (3 x 1) fitted initial distribution
%   P    (3 x 3 x nP) per-minute transition matrices, P(i,j,p) = P(j | i) in phase p
%   info struct: iterations, converged, startP (counts + 1, row-normalised)
%
%   The chain is the HMM with B fixed to the identity, fitted by HMMEMSTEP,
%   so missing minutes are marginalised exactly instead of breaking the
%   sequence. Start: transition counts + 1 per phase of the arrival minute
%   (pairs of adjacent observed minutes only). Stops when the relative
%   log-likelihood gain < Tol. On a complete record this is the count MLE.
%   A row with no expected transitions comes back as zeros (FILLEMPTYROWS
%   before simulating from it).
%
%   See also OBSERVABLEQUANTITIES, MODELENERGETICS, HMMEMSTEP.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    options.Tol (1,1) double {mustBePositive} = 1e-10
    options.MaxIter (1,1) double {mustBePositive, mustBeInteger} = 500
end
if numel(y) ~= numel(ph), error('sqw:length', 'y and ph must have the same length.'); end
I3 = eye(3);

t = (2:numel(y))';
a = y(t - 1); b = y(t);
use = a >= 1 & a <= 3 & b >= 1 & b <= 3;
P = ones(3, 3, nP) + accumarray([a(use), b(use), ph(t(use))], 1, [3 3 nP]);
P = P ./ sum(P, 2);
info.startP = P;
pi1 = ones(3, 1) / 3;
occ = ismember(1:3, y);
old = -Inf; info.converged = false;
for it = 1:options.MaxIter
    [pi1, P, B, llIn] = hmmEmStep(pi1, P, I3, y, ph);
    if any(abs(B(occ, :) - I3(occ, :)) > 1e-12, 'all')
        error('chainFitPhased:emissions', 'The identity emissions drifted: the chain fit is broken.');
    end
    if llIn - old < options.Tol * abs(llIn)
        info.converged = true;
        break
    end
    old = llIn;
end
info.iterations = it;
[~, ~, ll] = hmmForward(pi1, P, I3, y, ph);
end
