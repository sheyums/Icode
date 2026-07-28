function [ModelResults, fit_x] = shiftlognormal_MLE(x, mode)
% SHIFTLOGNORMAL_MLE Fit a constrained three-parameter shifted log-normal.
%
%   [ModelResults, fit_x] = shiftlognormal_MLE(x)
%   [ModelResults, fit_x] = shiftlognormal_MLE(x, mode)
%
%   fits the shifted log-normal model
%
%       X = shift + Y,
%       log(Y) ~ Normal(mu, sigma^2),
%
%   where X > shift.
%
%   IMPORTANT
%   ---------
%   The unrestricted three-parameter shifted log-normal likelihood is
%   unbounded as the shift approaches min(x) from below. This function
%   therefore constrains the shift to the data-derived interval
%
%       lowerShift = min(x) - range(x)
%
%       upperShift = min(x) - gapFraction*range(x)
%
%   with
%
%       gapFraction = 0.001.
%
%   These bounds ensure that the shift remains strictly less than min(x).
%   However, gapFraction is a regularization choice rather than a result
%   derived from maximum-likelihood theory. If the fitted shift is near a
%   bound, the function issues a warning.
%
%   INPUTS
%   ------
%   x
%       Real numeric vector containing the observations. The observations
%       may be negative. NaN and Inf values are removed before fitting.
%
%   mode
%       Optional character vector or string scalar.
%
%       If mode contains 'PLOT', the function plots a normalized histogram
%       and the fitted probability density.
%
%       'FIT' is accepted for compatibility, but the fitted PDF is always
%       calculated and returned in ModelResults.Fit.
%
%   OUTPUTS
%   -------
%   ModelResults
%       Structure containing the fitted model, parameter estimates,
%       likelihood, AIC, BIC, fit values, optimization diagnostics, and
%       shift-bound diagnostics.
%
%       Important fields include:
%
%           ModelResults.Params
%           ModelResults.Mu
%           ModelResults.Sigma
%           ModelResults.Shift
%           ModelResults.LogLike
%           ModelResults.AIC
%           ModelResults.BIC
%           ModelResults.FitX
%           ModelResults.Fit
%           ModelResults.Converged
%           ModelResults.ExitFlag
%           ModelResults.OptimizerOutput
%           ModelResults.Diagnostics
%
%   fit_x
%       Sorted x-coordinates corresponding to ModelResults.Fit.
%
%   EXAMPLE
%   -------
%       [results, fit_x] = shiftlognormal_MLE(x, 'PLOT');
%
%       plot(fit_x, results.Fit);
%
%   MODEL PARAMETERS
%   ----------------
%       Params(1) = mu
%       Params(2) = sigma
%       Params(3) = shift

    %% Constants and default settings

    gapFraction = 0.001;
    minimumSampleSize = 5;

    optimizerOptions = optimset( ...
        'Display', 'off', ...
        'TolX', 1e-8, ...
        'TolFun', 1e-8, ...
        'MaxIter', 500, ...
        'MaxFunEvals', 1000);

    %% Validate mode

    if nargin < 2 || isempty(mode)
        mode = '';
    end

    if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
        error('shiftlognormal_MLE:InvalidMode', ...
            'mode must be a character vector or string scalar.');
    end

    mode = upper(string(mode));

    %% Validate and clean the data

    validateattributes(x, {'numeric'}, ...
        {'real', 'vector', 'nonempty'}, mfilename, 'x', 1);

    x = x(:);

    finiteMask = isfinite(x);
    nRemoved = sum(~finiteMask);
    x = x(finiteMask);

    if nRemoved > 0
        warning('shiftlognormal_MLE:NonfiniteDataRemoved', ...
            '%d nonfinite observation(s) were removed.', nRemoved);
    end

    nSamples = numel(x);

    if nSamples < minimumSampleSize
        error('shiftlognormal_MLE:InsufficientData', ...
            ['At least %d finite observations are required, but only ' ...
             '%d were supplied.'], ...
            minimumSampleSize, nSamples);
    end

    min_x = min(x);
    max_x = max(x);
    dataRange = max_x - min_x;

    if ~isfinite(dataRange) || dataRange <= 0
        error('shiftlognormal_MLE:ConstantData', ...
            ['The finite observations must have a positive range. ' ...
             'A shifted log-normal distribution cannot be fitted to ' ...
             'constant data.']);
    end

    %% Construct the data-derived shift bounds

    minimumGap = gapFraction*dataRange;

    lowerShift = min_x - dataRange;
    upperShift = min_x - minimumGap;

    shiftBounds = [lowerShift, upperShift];

    % Confirm that the computed bounds are finite and properly ordered.
    if any(~isfinite(shiftBounds)) || lowerShift >= upperShift
        error('shiftlognormal_MLE:InvalidShiftBounds', ...
            'The data-derived shift bounds are invalid.');
    end

    % The upper bound must be strictly less than the smallest observation.
    if upperShift >= min_x
        error('shiftlognormal_MLE:InvalidUpperShift', ...
            ['The computed upper shift bound must be strictly less than ' ...
             'min(x).']);
    end

    fprintf('Fitting a constrained shifted log-normal model...\n');
    fprintf('  Shift bounds: [%.8g, %.8g]\n', ...
        lowerShift, upperShift);

    %% Optimize the profile likelihood over shift

    objectiveFunction = @(shift) ...
        neg_profile_loglike(x, shift);

    [bestShift, minNLL, exitFlag, optimizerOutput] = ...
        fminbnd( ...
            objectiveFunction, ...
            lowerShift, ...
            upperShift, ...
            optimizerOptions);

    converged = exitFlag > 0 && ...
                isfinite(bestShift) && ...
                isfinite(minNLL);

    if ~converged
        error('shiftlognormal_MLE:OptimizationFailed', ...
            ['The bounded shift optimization did not converge to a ' ...
             'finite solution. Exit flag: %d. Message: %s'], ...
            exitFlag, get_optimizer_message(optimizerOutput));
    end

    %% Recover mu and sigma at the fitted shift

    shiftedData = x - bestShift;

    if any(shiftedData <= 0) || any(~isfinite(shiftedData))
        error('shiftlognormal_MLE:InvalidShiftedData', ...
            ['The fitted shift produced nonpositive or nonfinite ' ...
             'shifted observations.']);
    end

    logData = log(shiftedData);

    bestMu = mean(logData);

    % The second argument of 1 uses normalization by N, which gives the
    % maximum-likelihood estimate of sigma for a normal distribution.
    bestSigma = std(logData, 1);

    if ~isfinite(bestMu) || ...
            ~isfinite(bestSigma) || bestSigma <= 0
        error('shiftlognormal_MLE:InvalidParameterEstimate', ...
            'The fitted mu or sigma estimate is invalid.');
    end

    %% Calculate likelihood and information criteria

    logLikelihood = -minNLL;

    % Mu, sigma, and shift are estimated from the data.
    nParameters = 3;

    % AIC = 2k - 2*log(L)
    AIC = 2*nParameters - 2*logLikelihood;

    % BIC = k*log(n) - 2*log(L)
    BIC = nParameters*log(nSamples) - 2*logLikelihood;

    %% Diagnose proximity to the shift bounds

    shiftIntervalWidth = upperShift - lowerShift;

    boundaryTolerance = max( ...
        10*optimizerOptions.TolX, ...
        1e-6*shiftIntervalWidth);

    distanceToLowerBound = bestShift - lowerShift;
    distanceToUpperBound = upperShift - bestShift;

    nearLowerBound = ...
        distanceToLowerBound <= boundaryTolerance;

    nearUpperBound = ...
        distanceToUpperBound <= boundaryTolerance;

    nearAnyBoundary = nearLowerBound || nearUpperBound;

    boundaryMessage = '';

    if nearLowerBound
        boundaryMessage = sprintf( ...
            ['The estimated shift is at or near the data-derived lower ' ...
             'bound of %.8g. The result may depend on how the lower ' ...
             'bound was defined.'], ...
            lowerShift);

        warning('shiftlognormal_MLE:ShiftNearLowerBound', ...
            '%s', boundaryMessage);

    elseif nearUpperBound
        boundaryMessage = sprintf( ...
            ['The estimated shift is at or near the data-derived upper ' ...
             'bound of %.8g. The result may depend strongly on the ' ...
             'selected gapFraction of %.8g.'], ...
            upperShift, gapFraction);

        warning('shiftlognormal_MLE:ShiftNearUpperBound', ...
            '%s', boundaryMessage);
    end

    %% Calculate fitted PDF and corresponding coordinates

    fit_x = sort(x);

    fitLogPDF = shifted_lognormal_logpdf( ...
        fit_x, bestMu, bestSigma, bestShift);

    fitPDF = exp(fitLogPDF);

    if any(~isfinite(fitPDF))
        warning('shiftlognormal_MLE:NonfiniteFittedDensity', ...
            'One or more fitted PDF values are nonfinite.');
    end

    %% Store results

    ModelResults = struct();

    ModelResults.Model = ...
        'Constrained shifted log-normal (mu, sigma, shift)';

    ModelResults.ParameterNames = ...
        {'mu', 'sigma', 'shift'};

    ModelResults.Params = ...
        [bestMu, bestSigma, bestShift];

    ModelResults.Mu = bestMu;
    ModelResults.Sigma = bestSigma;
    ModelResults.Shift = bestShift;

    ModelResults.LogLike = logLikelihood;
    ModelResults.AIC = AIC;
    ModelResults.BIC = BIC;

    ModelResults.FitX = fit_x;
    ModelResults.Fit = fitPDF;

    ModelResults.N = nSamples;
    ModelResults.NParameters = nParameters;

    ModelResults.Converged = converged;
    ModelResults.ExitFlag = exitFlag;
    ModelResults.OptimizerOutput = optimizerOutput;

    %% Store shift and boundary diagnostics

    ModelResults.Diagnostics = struct();

    ModelResults.Diagnostics.MinimumObservation = min_x;
    ModelResults.Diagnostics.MaximumObservation = max_x;
    ModelResults.Diagnostics.DataRange = dataRange;

    ModelResults.Diagnostics.GapFraction = gapFraction;
    ModelResults.Diagnostics.MinimumGap = minimumGap;

    ModelResults.Diagnostics.ShiftBounds = shiftBounds;
    ModelResults.Diagnostics.LowerShift = lowerShift;
    ModelResults.Diagnostics.UpperShift = upperShift;

    ModelResults.Diagnostics.ShiftToMinimumGap = ...
        min_x - bestShift;

    ModelResults.Diagnostics.BoundaryTolerance = ...
        boundaryTolerance;

    ModelResults.Diagnostics.DistanceToLowerBound = ...
        distanceToLowerBound;

    ModelResults.Diagnostics.DistanceToUpperBound = ...
        distanceToUpperBound;

    ModelResults.Diagnostics.NearLowerBound = ...
        nearLowerBound;

    ModelResults.Diagnostics.NearUpperBound = ...
        nearUpperBound;

    ModelResults.Diagnostics.NearAnyBoundary = ...
        nearAnyBoundary;

    ModelResults.Diagnostics.BoundaryMessage = ...
        boundaryMessage;

    ModelResults.Diagnostics.NonfiniteObservationsRemoved = ...
        nRemoved;

    %% Optional plot

    if contains(mode, "PLOT")
        figure('Color', 'w');

        histogram(x, ...
            'Normalization', 'pdf', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 1.5);

        hold on;

        plot(fit_x, fitPDF, ...
            'r-', ...
            'LineWidth', 2);

        xlabel('X value');
        ylabel('Probability density');

        title(sprintf( ...
            ['Shifted log-normal: \\mu = %.3g, ' ...
             '\\sigma = %.3g, shift = %.3g'], ...
            bestMu, bestSigma, bestShift));

        grid on;
        hold off;
    end

    %% Nested helper functions

    function nll = neg_profile_loglike(data, shift)
        % NEG_PROFILE_LOGLIKE Calculate the profile negative log-likelihood.
        %
        % For a specified shift, the maximum-likelihood estimates of mu
        % and sigma are calculated analytically. The resulting likelihood
        % is therefore a one-dimensional profile likelihood over shift.

        shiftedValues = data - shift;

        if any(shiftedValues <= 0) || ...
                any(~isfinite(shiftedValues))
            nll = Inf;
            return;
        end

        logShifted = log(shiftedValues);

        sigmaMLE = std(logShifted, 1);

        if ~isfinite(sigmaMLE) || sigmaMLE <= 0
            nll = Inf;
            return;
        end

        n = numel(data);

        % At the conditional MLEs of mu and sigma:
        %
        % NLL = n*log(sigma)
        %     + sum(log(x - shift))
        %     + n/2*(1 + log(2*pi)).
        nll = n*log(sigmaMLE) ...
            + sum(logShifted) ...
            + 0.5*n*(1 + log(2*pi));

        if ~isfinite(nll)
            nll = Inf;
        end
    end

    function logPDF = shifted_lognormal_logpdf( ...
            z, mu, sigma, shift)
        % SHIFTED_LOGNORMAL_LOGPDF Evaluate the shifted log-normal log-PDF.
        %
        % Values outside the support z > shift are assigned a log-density
        % of -Inf.

        shiftedZ = z - shift;

        logPDF = -Inf(size(z));

        valid = shiftedZ > 0 & isfinite(shiftedZ);

        if any(valid)
            logShiftedZ = log(shiftedZ(valid));

            logPDF(valid) = ...
                -logShiftedZ ...
                -log(sigma) ...
                -0.5*log(2*pi) ...
                -0.5*((logShiftedZ - mu)/sigma).^2;
        end
    end

    function message = get_optimizer_message(output)
        % GET_OPTIMIZER_MESSAGE Retrieve the termination message safely.

        if isstruct(output) && ...
                isfield(output, 'message') && ...
                ~isempty(output.message)

            message = output.message;
        else
            message = 'No optimizer message was returned.';
        end
    end
end
