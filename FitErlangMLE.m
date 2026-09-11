function H = FitErlangMLE(eventseries, xmin, varargin)
%FITERLANGMLE Left-truncated Erlang fit (gamma with INTEGER shape).
%
%   H = FITERLANGMLE(eventseries, xmin, SamplingInterval=dt, ...)
%   H = FITERLANGMLE(..., MaxShape=m)      default 12
%
%   The Erlang is the gamma restricted to integer shape. Because the shape
%   is discrete it cannot be optimized by fminsearch, so this sweeps
%   shape = 1..MaxShape, fits the scale at each, and keeps the best
%   log-likelihood. H.Params is [shape, scale] and H.ShapeSweep holds the
%   log-likelihood at every shape, so you can see whether the optimum is
%   interior or stuck at an endpoint.
%
%   TWO WARNINGS, BOTH STRUCTURAL.
%
%   The integer shape breaks the parameter count that AIC, AICc and BIC
%   assume. Those penalties are derived for continuous parameters; an
%   integer one is not worth a full degree of freedom. H.k is reported as 2
%   (shape and scale) which OVER-penalizes, so Erlang's AICc here is
%   conservative. Do not read a narrow Erlang defeat as decisive either way.
%
%   Erlang cannot represent a decreasing hazard AT ALL. Gamma's hazard
%   decreases only for shape < 1, and the smallest Erlang shape is 1, which
%   is the exponential (constant hazard). Bout durations generally have a
%   decreasing hazard, so this family is excluded from fitting them by
%   construction rather than by evidence. On a real 3415-bout sleep dataset
%   the free gamma wanted shape 0.134, and forcing shape 1 and 2 cost 362
%   and 1276 AICc units. Fit it for completeness of a sweep; do not expect
%   it to compete, and do not conclude anything biological from its defeat.
%
%   See also FITTRUNCATEDDISCRETEMLE, FITGAMMAMLE.

maxShape = 12;
keep = true(1, numel(varargin));
for ii = 1:numel(varargin)-1
    if (ischar(varargin{ii}) || isstring(varargin{ii})) && ...
            strcmpi(char(varargin{ii}), 'MaxShape')
        maxShape = varargin{ii+1};
        keep(ii:ii+1) = false;
    end
end
passThrough = varargin(keep);

best = []; sweep = nan(1, maxShape);
for k0 = 1:maxShape
    try
        Hk = FitTruncatedDiscreteMLE(eventseries, xmin, "gamma_fixedshape", ...
            passThrough{:}, 'Shape', k0);
        sweep(k0) = Hk.LogLik;
        if isempty(best) || (Hk.Success && Hk.LogLik > best.LogLik)
            best = Hk; best.ShapeInteger = k0;
        end
    catch err
        if strcmp(err.identifier, 'FitTruncatedDiscreteMLE:SamplingIntervalRequired') ...
                || strcmp(err.identifier, 'FitTruncatedDiscreteMLE:InvalidSamplingInterval')
            rethrow(err);
        end
    end
end
if isempty(best)
    error('FitErlangMLE:NoValidFit', ...
        'No integer shape in 1..%d produced a usable fit.', maxShape);
end

H = best;
H.Model = 'erlang';
H.ParamNames = {'shape','scale'};
H.Params = [best.ShapeInteger, best.Params(1)];
H.ParamSE = [NaN, best.ParamSE(1)];   % the integer shape has no standard error
H.k = 2;                              % see the header: this OVER-penalizes
H.AIC = 2*H.k - 2*H.LogLik;
if (H.n - H.k - 1) > 0
    H.AICc = H.AIC + (2*H.k*(H.k+1))/(H.n - H.k - 1);
else
    H.AICc = Inf;
end
H.BIC = H.k*log(H.n) - 2*H.LogLik;
H.ShapeSweep = sweep;
end
