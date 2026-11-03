function H = hmmPosteriorFluxEPR(y, ph, nP, fit, options)
%HMMPOSTERIORFLUXEPR  The HMM's counted cycle flux and EPR for one fly, averaged over posterior paths (sqw17.py's "hid"; Fig A1 triangles).
%
%   H = hmmPosteriorFluxEPR(y, ph, nP, fit)
%   H = hmmPosteriorFluxEPR(y, ph, nP, fit, NPaths=200, RandomSeed=0)
%
%   y, ph, nP  the fly's symbols and phases (SQWFLYINPUT)
%   fit        its LABELLED K = 3 20b fit (READSTOREDHMMFITS): pi, P, B, deg
%
%   H  (nP x 1) struct: hmm_deg (H1 degenerate, 20b), J_hmm = mean over NPaths
%      FFBS paths of JCYCLECOUNT (cycles/h), EPR_hmm = mean over paths of the
%      rule-of-three count EPR (EPR_COUNT_BASED 'bound', nats/h), and the model's
%      own J_hmm_model (MODELENERGETICS of P) and EPR_hmm_model (MODELENERGETICSR3: a one-way edge
%      filled by the rule of three from the MEAN path counts, r3_hmm = true; RuleOfThree=false: the
%      plain MODELENERGETICS, Inf kept, as before 2026-11-02). Counts on each path
%      are PAIRCOUNTS (phase of the arrival minute); hours = minutes in phase / 60.
%   EPR is averaged per path, E[f(N)], not computed on mean counts, as sqw2.decode.
%   Paths come from a private stream (HMMFFBS): compare with Python statistically.
%
%   See also HMMFFBS, PAIRCOUNTS, BUILDPERFLYDATA.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    fit (1,1) struct
    options.NPaths (1,1) double {mustBePositive, mustBeInteger} = 200
    options.RandomSeed (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    options.RuleOfThree (1,1) logical = true
end
Z = hmmFFBS(fit.pi(:), fit.P, fit.B, y, ph, NPaths=options.NPaths, RandomSeed=options.RandomSeed);
hours = accumarray(ph, 1, [nP 1]) / 60;
J = zeros(options.NPaths, nP); E = J; Nm = zeros(3, 3, nP);
for k = 1:options.NPaths
    N = pairCounts(Z(k, :)', ph, nP); Nm = Nm + N / options.NPaths;
    J(k, :) = jCycleCount(N, hours)';
    for p = 1:nP, E(k, p) = epr_count_based(N(:, :, p), hours(p) * 3600, 'bound'); end
end
H = repmat(struct('hmm_deg', logical(fit.deg), 'J_hmm', NaN, 'EPR_hmm', NaN, 'J_hmm_model', NaN, ...
    'EPR_hmm_model', NaN, 'r3_hmm', false), nP, 1);
for p = 1:nP
    e = modelEnergetics(fit.P(:, :, p));
    H(p).J_hmm = mean(J(:, p)); H(p).EPR_hmm = mean(E(:, p));
    H(p).J_hmm_model = e.J; H(p).EPR_hmm_model = e.EPR;
    if options.RuleOfThree        % one-way edges filled from the mean posterior-path counts (sqw_r3.py)
        er = modelEnergeticsR3(fit.P(:, :, p), Nm(:, :, p)); H(p).EPR_hmm_model = er.EPR; H(p).r3_hmm = er.r3;
    end
end
end
