function H = FitHyperErlangMLE(eventseries, xmin, varargin)
%FITHYPERERLANGMLE Left-truncated hyper-Erlang: sum_j q_j Erlang(m_j, lam_j).
%
%   H = FITHYPERERLANGMLE(eventseries, xmin, SamplingInterval=dt, ...)
%   H = FITHYPERERLANGMLE(..., Components=3, MaxShape=6)
%   H = FITHYPERERLANGMLE(..., Shapes=[1 1 2])     one fixed shape vector
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
keep = true(1, numel(varargin));
for ii = 1:numel(varargin)-1
    if ~(ischar(varargin{ii}) || isstring(varargin{ii})), continue; end
    switch lower(char(varargin{ii}))
        case 'components', nComp = varargin{ii+1};       keep(ii:ii+1) = false;
        case 'maxshape',   maxShape = varargin{ii+1};    keep(ii:ii+1) = false;
        case 'shapes',     fixedShapes = varargin{ii+1}; keep(ii:ii+1) = false;
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
for m = 1:maxShape
    shp = [ones(1, nComp-1), m];
    try
        Hm = FitTruncatedDiscreteMLE(eventseries, xmin, "hyper_erlang", ...
            passThrough{:}, 'Shapes', shp);
        sweep(m) = Hm.LogLik;
        % A shape whose fit the guard rejects must not win the sweep: its
        % likelihood is real but its parameters are not identified, so
        % picking it would report estimates that do not mean anything.
        usable = Hm.Success && Hm.Diagnostics.GuardOK;
        if usable && (isempty(best) || Hm.LogLik > best.LogLik)
            best = Hm; bestM = m;
        end
    catch err
        if strcmp(err.identifier, 'FitTruncatedDiscreteMLE:SamplingIntervalRequired') ...
                || strcmp(err.identifier, 'FitTruncatedDiscreteMLE:InvalidSamplingInterval')
            rethrow(err);
        end
    end
end
if isempty(best)
    error('FitHyperErlangMLE:NoValidFit', ...
        ['No shape vector in 1..%d produced an identified fit at %d ' ...
         'components. Try fewer components, or inspect a single Shapes ' ...
         'vector directly.'], maxShape, nComp);
end

H = best;
H.Model = 'hyper_erlang';
H.Shapes = [ones(1, nComp-1), bestM];
H.ShapeSweep = sweep;
H.SelectedShape = bestM;
end
