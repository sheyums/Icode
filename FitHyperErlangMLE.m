function H = FitHyperErlangMLE(eventseries, xmin, varargin)
%FITHYPERERLANGMLE Left-truncated hyper-Erlang: sum_j q_j Erlang(m_j, lam_j).
%
%   H = FITHYPERERLANGMLE(eventseries, xmin, SamplingInterval=dt, ...)
%   H = FITHYPERERLANGMLE(..., Components=3, MaxShape=6)
%   H = FITHYPERERLANGMLE(..., Shapes=[1 1 2])     one fixed shape vector
%   H = FITHYPERERLANGMLE(..., WarmStart=true, SweepStarts=3)
%   H = FITHYPERERLANGMLE(..., SeedRates=r, SeedWeights=q)
%
%   THE EXPONENTIAL-NATIVE WAY TO GET A NON-MONOTONE HAZARD. A
%   hyperexponential arranges its phases in PARALLEL -- enter one of K
%   states, leave at that state's own rate -- and that topology forces a
%   completely monotone hazard at EVERY K. No number of components will
%   produce a hump. An Erlang arranges phases in SERIES, must pass through
%   m stages before leaving, and its hazard RISES. Allow both in one
%   mixture and the hazard can fall, rise and fall again, while every
%   parameter remains a rate or a branch probability.
%
%   IT NESTS THE HYPEREXPONENTIAL EXACTLY. Shapes = ones(1,K) reproduces
%   FITHYPEREXPONENTIALMLE at order K to the digit, which makes "are
%   memoryless states enough?" a constraint on this one model rather than a
%   comparison between two. The sweep always includes that vector, so the
%   profile over m tells you directly what the extra stages bought.
%
%   INTERPRETATION. The fit is a phase-type of order sum(m_j), so it keeps
%   a Markov reading: a branch with m=2 says a bout in that branch passes
%   through two sequential sub-stages before it can end -- a refractory or
%   cumulative process rather than a memoryless one. Hyper-Erlangs are
%   dense in the distributions on [0,inf) (Tijms 1994, Asmussen 2003), so a
%   hyper-Erlang that still will not fit is evidence of something
%   structural rather than of a family too small.
%
%   An order-2 phase-type cannot hump at all: PH(2) hazards are monotone.
%   So the smallest shape vector that can produce one has sum(m_j) >= 3.
%
%   THE SWEEP. Integer shapes cannot be optimized by fminsearch, so this
%   sweeps shape vectors [1 ... 1 m] for m = 1..MaxShape over Components
%   branches -- every branch memoryless except one that needs m stages.
%   m=1 is the plain hyperexponential. Pass Shapes to fix one vector and
%   skip the sweep. H.ShapeSweep holds the log-likelihood at every m, so
%   you can see whether the optimum is interior or stuck at an endpoint.
%
%   COST, AND THE WARM START. The sweep is MaxShape independent fits, each
%   a full multistart, so at the engine's default it is MaxShape x 24
%   optimizations of a 2J-1 parameter mixture to produce ONE row --
%   measured at 379 s for m <= 4 against 49 s for seven other families
%   combined. WarmStart=true (the default) chains them instead: the
%   optimum at m-1 is carried into m with the series branch's rate
%   rescaled by m/(m-1), which holds that branch's MEAN fixed, and only
%   SweepStarts (default 3) cold starts are run alongside. The cold starts
%   are reduced, never removed, so the warm point can only add a candidate
%   optimum. Set WarmStart=false for the original independent sweep.
%
%   SeedRates/SeedWeights supply an external starting point for m=1 --
%   COMPAREBOUTMODELS passes the hyperexponential fit it has already
%   computed at the same order, which is the m=1 member of this very
%   family and therefore free.
%
%   READ THE SWEEP, NOT ONLY THE AICc. The integer shape breaks the
%   parameter count the information criteria assume: their penalties are
%   derived for continuous parameters, and an integer one is not worth a
%   full degree of freedom. H.k counts 2J-1 -- the rates and the free
%   weights -- and charges the shape NOTHING, which UNDER-penalizes, the
%   opposite of the convention FITERLANGMLE takes. Neither is right. A
%   shape chosen by maximizing the likelihood over a sweep has cost
%   something, and no information criterion here accounts for it, so treat
%   a narrow hyper-Erlang victory as undecided and look at H.ShapeSweep for
%   whether the likelihood really preferred m > 1 or merely drifted.
%
%   Check H.Diagnostics.GuardOK before reading the estimates: it is false
%   when a component holds fewer than MinMixtureCount observations, or when
%   two components sharing a shape have converged on the same rate. Both
%   leave parameters unidentified.
%
%   See also FITHYPEREXPONENTIALMLE, FITTRUNCATEDDISCRETEMLE,
%   FITWEIBULLMIXTUREMLE, FITERLANGMLE.

nComp = 3; maxShape = 6; fixedShapes = [];
warmStart = true; sweepStarts = 3; seedRates = []; seedWeights = [];
keep = true(1, numel(varargin));
for ii = 1:numel(varargin)-1
    if ~(ischar(varargin{ii}) || isstring(varargin{ii})), continue; end
    switch lower(char(varargin{ii}))
        case 'components', nComp = varargin{ii+1};       keep(ii:ii+1) = false;
        case 'maxshape',   maxShape = varargin{ii+1};    keep(ii:ii+1) = false;
        case 'shapes',     fixedShapes = varargin{ii+1}; keep(ii:ii+1) = false;
        case 'warmstart',  warmStart = varargin{ii+1};   keep(ii:ii+1) = false;
        case 'sweepstarts',sweepStarts = varargin{ii+1}; keep(ii:ii+1) = false;
        case 'seedrates',  seedRates = varargin{ii+1};   keep(ii:ii+1) = false;
        case 'seedweights',seedWeights = varargin{ii+1}; keep(ii:ii+1) = false;
    end
end
passThrough = varargin(keep);

if ~isempty(fixedShapes)
    H = FitTruncatedDiscreteMLE(eventseries, xmin, "hyper_erlang", ...
        passThrough{:}, 'Shapes', fixedShapes);
    H.Model = 'hyper_erlang';
    H.Shapes = sort(fixedShapes(:).');
    H.ShapeSweep = [];
    return
end

best = []; bestM = NaN; sweep = nan(1, maxShape);
bestAny = []; bestAnyM = NaN;   % best fit regardless of the guard
prevZ = [];                      % the previous shape's optimum, in z
for m = 1:maxShape
    shp = [ones(1, nComp-1), m];
    try
        % WARM START. The fit at m stages is a short step from the fit at
        % m-1 provided the branch MEAN is held fixed, because the mean is
        % what the data determine: an Erlang(m, lam) branch has mean
        % m/lam, so carrying lam forward unchanged would make the branch
        % suddenly m/(m-1) times slower and throw away the step. Rescaling
        % lam -> lam*m/(m-1) on the series branch moves the shape while
        % leaving the fitted timescale where the likelihood put it.
        %
        % That the mean is the identified quantity and the stage count is
        % not is not an assumption here -- it is what the engine's own
        % test 26 measures, where data generated with m=3 are fitted at
        % m=2 with the survival still within 0.015 of truth.
        %
        % The cold starts are REDUCED, not removed (sweepStarts of them),
        % so the warm point is an addition to the search, not a
        % substitute: an optimum somewhere else can still be found.
        warm = {};
        if warmStart
            if m == 1 && ~isempty(seedRates) && numel(seedRates) == nComp
                % SeedWeights are OBSERVED weights, the convention every
                % mixture here reports in. The model's own weights are
                % untruncated, so back-transform first:
                %   w_j  proportional to  q_j / S_j(xmin)
                % and at m=1 every component is exponential, so
                % S_j(xmin) = exp(-rate_j*xmin) and the gain is
                % exp(+rate_j*xmin). That factor reaches 240x in these
                % data, so it is taken in LOG space and re-centred before
                % exponentiating; done naively it overflows to Inf and the
                % seed becomes a wall instead of a hint.
                lw = log(max(seedWeights(:).', realmin)) + seedRates(:).'*xmin;
                lw = lw - max(lw);
                wSeed = exp(lw); wSeed = wSeed / sum(wSeed);
                z0 = heZ(seedRates, wSeed, nComp);
                if ~isempty(z0), warm = {'StartZ', z0, 'nStarts', sweepStarts}; end
            elseif m > 1 && ~isempty(prevZ)
                z0 = prevZ;
                z0(nComp) = z0(nComp) + log(m / (m-1));   % rates are logged
                warm = {'StartZ', z0, 'nStarts', sweepStarts};
            end
        end
        Hm = FitTruncatedDiscreteMLE(eventseries, xmin, "hyper_erlang", ...
            passThrough{:}, warm{:}, 'Shapes', shp);
        sweep(m) = Hm.LogLik;
        % Chain from this shape's optimum whether or not its GUARD passed:
        % a guard-rejected fit is still the likelihood's maximum at that
        % shape, and so still the best place to start the next one. Only
        % the WINNER is filtered on the guard, below.
        if Hm.Success
            % Theta, NOT Params: the reported weights are the OBSERVED
            % ones, and feeding those back as if they were the mixing
            % weights would start the next shape at a point the model
            % never occupied.
            zc = heZ(Hm.Theta(1:nComp), Hm.Theta(nComp+1:2*nComp), nComp);
            if ~isempty(zc), prevZ = zc; end
        end
        % A shape whose fit the guard rejects must not win the sweep: its
        % likelihood is real but its parameters are not identified, so
        % picking it would report estimates that do not mean anything.
        usable = Hm.Success && Hm.Diagnostics.GuardOK;
        if usable && (isempty(best) || Hm.LogLik > best.LogLik)
            best = Hm; bestM = m;
        end
        if Hm.Success && (isempty(bestAny) || Hm.LogLik > bestAny.LogLik)
            bestAny = Hm; bestAnyM = m;      % guard verdict ignored here
        end
    catch err
        if strcmp(err.identifier, 'FitTruncatedDiscreteMLE:SamplingIntervalRequired') ...
                || strcmp(err.identifier, 'FitTruncatedDiscreteMLE:InvalidSamplingInterval')
            rethrow(err);
        end
    end
end
if isempty(best)
    % Every shape was rejected by the guard. That is a real answer, not a
    % missing one -- typically Components is larger than the data support,
    % so whichever shape is tried one component ends up holding almost
    % none of the retained bouts. Return the best of them with its guard
    % verdict INTACT (GuardOK=false, with the reason) rather than
    % erroring, so a comparison table shows the row and says why it cannot
    % be used, which is how every other family here behaves. A caller that
    % ranks on the guard, as COMPAREBOUTMODELS does, will not let it win.
    if isempty(bestAny)
        error('FitHyperErlangMLE:NoValidFit', ...
            ['No shape vector in 1..%d could be fitted at all at %d ' ...
             'components.'], maxShape, nComp);
    end
    best = bestAny; bestM = bestAnyM;
end

H = best;
H.Model = 'hyper_erlang';
H.Shapes = [ones(1, nComp-1), bestM];
H.ShapeSweep = sweep;
H.SelectedShape = bestM;
end

% ------------------------------------------------------------------------
function z = heZ(rates, q, J)
%HEZ Natural parameters -> the engine's internal coordinates for
%"hyper_erlang": z = [log(rates), logits], the last logit pinned at 0 by
%the softmax, exactly as heUnpack reads them. Returns empty if the point
%cannot be expressed -- a vanished component has no finite logit, and a
%start built from one would be a wall, not a hint.
z = [];
if numel(rates) ~= J || numel(q) ~= J, return; end
rates = rates(:).'; q = q(:).';
if ~all(isfinite(rates)) || ~all(rates > 0), return; end
if ~all(isfinite(q)) || any(q <= 0), return; end
q = q / sum(q);
if J == 1
    z = log(rates);
else
    z = [log(rates), log(q(1:J-1)) - log(q(J))];
end
if ~all(isfinite(z)), z = []; end
end
