function [EPR, J_cycle, A_eff, nEdge] = epr_count_based(N, nSteps, rule, PMIN)
%EPR_COUNT_BASED  Schnakenberg entropy production from global transition counts.
%
% SYNTAX
%   [EPR, J_cycle, A_eff, nEdge] = epr_count_based(N, nSteps)          % 'bound'
%   [EPR, J_cycle, A_eff, nEdge] = epr_count_based(N, nSteps, rule, PMIN)
%
%   N        3 x 3 directed transition counts on the compressed path. Either
%            from COMPRESSED_PAIR_COUNTS, or assembled from a
%            COMPUTE_HMM_ENERGETICS fly struct as
%               N = [0 f.N12 f.N13; f.N21 0 f.N23; f.N31 f.N32 0]
%   nSteps   total valid time steps (= f.total_hours*3600 for that struct)
%   rule     how to treat an edge pair with a ZERO directed count:
%              'bound'  (DEFAULT) cap the unobserved rate at 3/T, the 95%
%                       rule-of-three limit, giving a LOWER BOUND on EPR
%              'drop'   or 0 -- the published convention: the pair contributes
%                       nothing, i.e. +Inf is replaced by 0
%              numeric a > 0 -- shrinkage, J = (N + a)/T
%   PMIN     observability floor, default 1e-7, used only by 'drop'
%
%   EPR      nats/hr      J_cycle  transitions/hr (INVARIANT to rule)
%   A_eff    EPR/|J_cycle|, nats/cycle        nEdge  pairs that contributed
%
% ------------------------------------------------------------------------
% WHY THE DEFAULT IS A BOUND
%
%   The Schnakenberg term for one edge is (Jf - Jr)*log(Jf/Jr). As Jr -> 0 this
%   DIVERGES: a perfectly irreversible step dissipates infinitely per unit flux.
%   A zero reverse count is therefore not a missing measurement, it is the
%   statement "the reverse rate is below my detection limit", and the correct
%   treatment is a bound, not a substitution.
%
%   P(observe 0 | rate lambda over T) = exp(-lambda*T), so requiring that to be
%   at least 0.05 bounds the unobserved rate by lambda < 3/T (this is the 95%
%   one-sided Poisson limit, ln(20) = 2.996; the "rule of three" of the medical
%   statistics literature).
%
%   THE CAP MUST ENTER BOTH THE PREFACTOR AND THE LOG. The term as a function
%   of the unobserved rate, f(Jr) = (k/T - Jr)*log((k/T)/Jr), diverges as
%   Jr -> 0 and decreases monotonically to zero at Jr = k/T. So the minimum
%   over the permitted interval sits at Jr = min(3/T, k/T):
%
%       k <= 3   ->  0, the counts are consistent with detailed balance
%       k >  3   ->  ((k-3)/T) * log(k/3)
%
%   An earlier version applied the cap only inside the logarithm while leaving
%   the raw zero in the prefactor. That flips the sign of the term whenever
%   k <= 3 and produced a NEGATIVE EPR for one sss fly, which is impossible --
%   every Schnakenberg term is a product of two factors of the same sign. The
%   bug survived six rounds of verification because every check was a
%   consistency check ('drop' still reproduces the published values', 'two code
%   paths agree') and none asserted the physical constraint EPR >= 0. Worth
%   remembering: check invariants, not just agreement.
%
%   This treatment preserves the Methods' framing of EPR_obs as a lower bound,
%   now in a second and stronger sense than the data-processing inequality.
%
%   The published 'drop' rule instead replaces the divergence by ZERO -- the
%   maximum possible value by the minimum. Where that bites is not uniform:
%   share of edge pairs with one directed count zero, measured 2026-10-15 --
%
%       wt 0%,  per0 DD 0%,  per0 LD 2%,  fmn 53%,  sss 51%
%
%   and in fmn/sss those zeros are the REVERSE LEGS OF A STRONGLY DRIVEN CYCLE.
%   Directed traffic in fmn: H2->H1 1310, H1->H3 1321, H3->H2 2152 (52-53 of 53
%   flies) against reverse legs H3->H1 97 and H1->H2 87 (9-10 of 53 flies).
%   Those flies run the Sleep->Wake->Reg->Sleep cycle almost one-way, which is
%   exactly the condition for LARGE dissipation -- so dropping the edges scores
%   the most irreversible flies in the dataset as the least dissipative.
%
% WHAT THE RULE CHANGES  (the paper's five genotypes with the EXTRA-DAYS fmn
%                         and sss records; reproduce with COMPARE_EPR_ESTIMATORS)
%
%   median EPR   drop   BOUND         |  L          drop   BOUND   P(L>L_wt) bnd
%     wt        0.579   0.579         |  wt        0.228   0.228
%     per0-pos  0.437   0.437         |  per0-pos  0.256   0.256      0.870
%     per0-neg  2.524   2.703         |  per0-neg  0.257   0.185      0.081
%     fmn       0.725   3.493         |  fmn       0.217   0.111      0.000
%     sss       1.342   2.599         |  sss       0.199   0.126      0.000
%
%   On the SHORT fmn/sss records the same table read fmn 0.405 -> 1.533 with
%   L 0.401 (drop) and 0.148 (bound), sss 0.348 -> 1.107 with L 0.474 and 0.173.
%   The published claim "fmn and sss show significantly HIGHER L than wt" held
%   only there, and only under 'drop': on the longer records 'drop' itself gives
%   P = 0.386 and 0.190, i.e. no difference. The BOUND gives the same verdict at
%   both record lengths (significantly lower, P <= 0.001), which is the strongest
%   argument for it -- its answer does not depend on recording duration.
%
%   wt and per0-pos are IDENTICAL under all four rules, because they have no
%   zero counts at all -- the rule acts only where it should, which is the
%   built-in check that it is doing what it claims. 'bound' is the most
%   conservative of the three corrected rules in every row.
%
%   Bootstrap P(L > L_wt). The published claim is "fmn and sss show
%   significantly HIGHER L than wt"; only 'drop' gives that.
%
%                  drop   BOUND  shr.5   shr1
%       per0-pos   0.860  0.866  0.875  0.885    ns under every rule
%       per0-neg   0.653  0.086  0.075  0.089    ns, but marginal once corrected
%       fmn        0.997  0.000  0.000  0.000    REVERSES
%       sss        0.999  0.001  0.000  0.000    REVERSES
%
%   Under 'bound': fmn L = 0.148 [0.134 0.157] and sss 0.173 [0.141 0.195]
%   against wt 0.228 [0.206 0.249] -- non-overlapping. A_eff also reverses, fmn
%   3.80 and sss 3.00 against wt 1.65, where 'drop' put them at 0.84 and 0.79.
%
%   Note per0-neg: its EPR moves only 7% but its L moves 28% and its P goes
%   0.653 -> 0.086. The verdict (no difference from wt) holds under all four
%   rules, but it is marginal rather than comfortable. That is the
%   leverage sensitivity of the origin-forced fit on x = J^2, not the estimator.
%
%   'bound' is also LOWER-VARIANCE than 'drop': its fmn interval is several
%   times narrower. Under 'drop' each fly contributes a different number of
%   edges, so the scatter partly measures the gate rather than the fly.
%
% HOW TO REPORT IT
%   For fmn and sss the defensible statement is not "EPR is not estimable" but
%   "EPR is at least 1.53 and 1.11 nats/hr and could be considerably higher,
%   because the reverse transitions are too rare to measure."
%
% CAVEAT
%   'bound' has no free parameter, which is its main advantage over shrinkage
%   (where fmn EPR swings 15% between a = 0.5 and a = 1). It is conservative by
%   construction: a conclusion that survives at the bound survives under any
%   less conservative treatment.
%
% See also COMPARE_EPR_ESTIMATORS, COMPRESSED_PAIR_COUNTS,
%          EXPORT_FIG6_THERMODYNAMICS_FORORIGIN, COMPUTE_HMM_ENERGETICS

if nargin < 3 || isempty(rule),  rule = 'bound'; end
if nargin < 4 || isempty(PMIN),  PMIN = 1e-7;    end

EPR = NaN; J_cycle = NaN; A_eff = NaN; nEdge = 0;
if ~isfinite(nSteps) || nSteps <= 0, return; end

Th = nSteps / 3600;
J_cycle = (N(3,2) - N(2,3)) / Th;     % pseudocounts and caps cancel here

isDrop = (ischar(rule) || isstring(rule)) && strcmpi(rule,'drop');
isBound = (ischar(rule) || isstring(rule)) && strcmpi(rule,'bound');
if isnumeric(rule) && rule == 0, isDrop = true; end

pairs = [1 2; 1 3; 2 3];
e = 0; c = 0;
for q = 1:3
    i = pairs(q,1); j = pairs(q,2);
    u = N(i,j); v = N(j,i);
    if isDrop
        if u/nSteps <= PMIN || v/nSteps <= PMIN, continue; end
        e = e + (u/Th - v/Th)*log(u/v); c = c + 1;
    elseif isBound
        if u == 0 && v == 0, continue; end      % no evidence either way
        if u > 0 && v > 0
            e = e + (u/Th - v/Th)*log(u/v);
        else
            % One direction unobserved. Its rate lies in [0, 3/T] by the rule
            % of three, and the term (k/T - Jr)*log((k/T)/Jr) DECREASES
            % monotonically in Jr, reaching 0 at Jr = k/T. So the minimum over
            % the allowed interval sits at Jr = min(3/T, k/T):
            %   k <= 3  ->  0, the counts are consistent with reversibility
            %   k >  3  ->  ((k-3)/T)*log(k/3)
            % The cap must enter BOTH the prefactor and the log; applying it to
            % only one of them flips the sign of the term for small k and can
            % make EPR negative, which is impossible.
            k = max(u, v);
            if k > 3
                e = e + ((k - 3)/Th) * log(k/3);
            end
        end
        c = c + 1;
    else
        a = double(rule);
        e = e + ((u+a)/Th - (v+a)/Th)*log((u+a)/(v+a)); c = c + 1;
    end
end
nEdge = c;
if c > 0, EPR = e; end
A_eff = EPR / max(abs(J_cycle), realmin);
end
