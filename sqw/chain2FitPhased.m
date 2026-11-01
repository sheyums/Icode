function [ll, pi1, P, info] = chain2FitPhased(y, ph, nP, options)
%CHAIN2FITPHASED  Second-order Markov chain on the observed S/Q/W symbols, per phase (sqw_ref.mc2_fit; PRESPEC 20o).
%
%   [ll, pi1, P, info] = chain2FitPhased(y, ph, nP)
%   [ll, pi1, P, info] = chain2FitPhased(y, ph, nP, Tol=1e-10, MaxIter=1000)
%
%   y    (T x 1) symbols 1 = Q, 2 = W, 3 = S, 0 = missing, 4 = unknown wake
%   ph   (T x 1) phase of each minute, 1..nP (arrival-minute convention)
%
%   ll   log-likelihood of the fitted chain
%   pi1  (9 x 1) fitted initial distribution over pair states
%   P    (9 x 9 x nP) transition matrices on pair states k = 3 (y_{t-1} - 1) + y_t;
%        P(k,l,p) is nonzero only where Allow(k,l) (PAIRCHAINSTRUCTURE), so
%        P(y_t | y_{t-2}, y_{t-1}) = P(3 (y_{t-2} - 1) + y_{t-1}, 3 (y_{t-1} - 1) + y_t, p)
%   info struct: iterations, converged, startP, B (the fixed 9 x 3 emissions),
%        Allow
%
%   P(y_t | y_{t-2}, y_{t-1}) as a 9-state chain on pairs with FIXED emissions
%   (pair k emits its second symbol) and structural zeros, fitted by HMMEMSTEP
%   so missing minutes are marginalised exactly. Start: triple counts + 1 on the
%   allowed transitions, per phase of the arrival minute (triples of observed
%   minutes only), pi uniform. Stops when the relative log-likelihood gain
%   < Tol. A row with no expected transitions comes back as zeros
%   (FILLEMPTYROWS(P, Allow) before simulating).
%
%   See also CHAINFITPHASED, PAIRCHAINSTRUCTURE, HMMEMSTEP.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    options.Tol (1,1) double {mustBePositive} = 1e-10
    options.MaxIter (1,1) double {mustBePositive, mustBeInteger} = 1000
end
if numel(y) ~= numel(ph), error('sqw:length', 'y and ph must have the same length.'); end
[B2, Allow] = pairChainStructure();

t = (3:numel(y))';
a = y(t - 2); b = y(t - 1); c = y(t);
use = a >= 1 & a <= 3 & b >= 1 & b <= 3 & c >= 1 & c <= 3;
k = 3 * (a(use) - 1) + b(use); l = 3 * (b(use) - 1) + c(use);
P = repmat(double(Allow), 1, 1, nP) + accumarray([k, l, ph(t(use))], 1, [9 9 nP]);
P = P ./ sum(P, 2);
info.startP = P;
pi1 = ones(9, 1) / 9;
old = -Inf; info.converged = false;
for it = 1:options.MaxIter
    [pi1, P, ~, llIn] = hmmEmStep(pi1, P, B2, y, ph);
    if any(P(~repmat(Allow, 1, 1, nP)) ~= 0)
        error('chain2FitPhased:structure', 'A structural zero became nonzero: the pair-chain fit is broken.');
    end
    if llIn - old < options.Tol * abs(llIn)
        info.converged = true;
        break
    end
    old = llIn;
end
info.iterations = it; info.B = B2; info.Allow = Allow;
[~, ~, ll] = hmmForward(pi1, P, B2, y, ph);
end
