function ModelResults = shiftlognormal_MLE(x, mode)
% SHIFTLOGNORMAL_MLE Fit a constrained three-parameter shifted log-normal.
%
%   ModelResults = shiftlognormal_MLE(x)
%   ModelResults = shiftlognormal_MLE(x, mode)
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
%       upperShift = min(x) - gapFraction*median(x - min(x))
%
%   with gapFraction = 0.001 (range(x) replaces the median when more than
%   half the observations tie at min(x)).
%
%   Upper bound: gapFraction controls the minimum permitted gap between
%   the fitted shift and the smallest observation. This is a regularization
%   choice, not a result derived from maximum-likelihood theory. If the
%   fitted shift approaches this bound, the estimate is sensitive to
%   gapFraction and a warning is issued. The gap is scaled by the median
%   distance above min(x), not by the range: the range is dominated by
%   the right tail, and for heavy tails 0.001*range(x) exceeded the true
%   distance from the shift to min(x), putting the true shift OUTSIDE the
%   allowed interval (measured Nov 2026, n = 2000, true shift 3: sigma = 1
%   clipped 10/20 fits, sigma = 1.5 clipped 20/20 with a median error of
%   0.45; with the median scale, 0/20 and 0/20, median error 0.006).
%
%   The search runs over t = log(min(x) - shift), and the shifted data are
%   formed as (x - min(x)) + gap. The fit is therefore equivariant under
%   x -> a*x + b: shift -> a*shift + b, mu -> mu + log(a), sigma unchanged.
%
%   Lower bound: the shift cannot be further than one full data range below
%   the minimum observation. This is a heuristic. If the true shift lies
%   further left, the optimizer will clip to the lower boundary and a
%   warning is issued. In that case, widen the lower bound manually or
%   interpret the result with caution.
%
%   If the fitted shift is near either bound, the result should be
%   interpreted with caution and ModelResults.Diagnostics inspected.
%
%   INPUTS
%   ------
%   x
%       Real numeric vector containing the observations. The observations
%       may be negative. NaN and Inf values are removed before fitting.
%
%   mode
%       Optional character vector or string scalar. Multiple keywords may
%       be combined in a single string (e.g., 'PLOT VERBOSE').
%
%       'PLOT'    - Plot a normalized histogram with the fitted PDF overlay.
%       'VERBOSE' - Print progress and bound information to the command
%                   window. 'PLOT' implies 'VERBOSE'.
%       'FIT'     - Accepted for compatibility; has no additional effect.
%
%   OUTPUTS
%   -------
%   ModelResults
%       Structure with fields:
%           .Model             Description string
%           .ParameterNames    {'mu', 'sigma', 'shift'}
%           .Params            [mu, sigma, shift]
%           .Mu                Fitted mu
%           .Sigma             Fitted sigma (MLE, normalized by n)
%           .Shift             Fitted shift
%           .LogLike           Profile log-likelihood at the MLE
%           .AIC               Akaike information criterion (2k - 2 logL)
%           .BIC               Bayesian information criterion
%                              (k*log(n) - 2 logL)
%           .FitX              Sorted x values used to evaluate the PDF
%           .Fit               Fitted PDF evaluated at FitX
%           .N                 Number of observations used after cleaning
%           .NParameters       Number of free parameters (3)
%           .Converged         True if the refinement that produced the
%                              returned solution converged. Converged
%                              candidates are preferred over nonconverged
%                              candidates. If no refinement converged, the
%                              finite nonconverged candidate with the
%                              lowest NLL is returned and Converged is
%                              false.
%           .ExitFlag          fminbnd exit flag for the refinement that
%                              produced the returned solution
%           .OptimizerOutput   fminbnd output struct for the refinement
%                              that produced the returned solution
%           .Diagnostics       Struct with bound, grid, and data diagnostics.
%                              NearLowerBound / NearUpperBound are judged
%                              in t = log(gap), and BoundaryTolerance is in
%                              t (1e-6 of the t interval). GapScale is the
%                              median(x - min(x)) the upper bound uses.
%
%   WARNINGS
%   --------
%       shiftlognormal_MLE:SmallSample          n < 15 (see LIMITATIONS)
%       shiftlognormal_MLE:ShiftNearUpperBound  shift on the min(x) spike
%       shiftlognormal_MLE:ShiftNearLowerBound  shift clipped on the left
%       shiftlognormal_MLE:PartialConvergence   no refinement converged
%       shiftlognormal_MLE:NonfiniteDataRemoved NaN/Inf dropped
%
%   LIMITATIONS (measured Nov 2026)
%   -----------
%   - The likelihood is unbounded as shift -> min(x), so at small n the
%     spike at the upper bound often beats any interior optimum: with
%     log-normal data (sigma = 0.5) it won in 134/200 fits at n = 5,
%     34/200 at n = 8, 9/200 at n = 10, 5/200 at n = 12 and 0/200 from
%     n = 15. Hence the SmallSample warning below 15.
%   - Data with no right skew (Normal, uniform, left-skewed) have no
%     interior optimum: the profile keeps rising toward the Normal limit,
%     shift -> -Inf, and the fit clips to the lower bound with a warning.
%     Read such a fit as "no evidence of a shifted log-normal shape".
%   - No standard errors or confidence intervals are returned.
%
%   EXAMPLE
%   -------
%       results = shiftlognormal_MLE(x, 'PLOT VERBOSE');
%       plot(results.FitX, results.Fit);
%
%   MODEL PARAMETERS
%   ----------------
%       Params(1) = mu
%       Params(2) = sigma
%       Params(3) = shift

%% Original written by Sheyum. Hardened by ChatGPT and Claude.
%  Audited against an independent profile likelihood (Claude Code,
%  Nov 2026): the search moved to t = log(min(x) - shift), fixing a loss
%  of equivariance under large offsets (x + 1e8 moved the shift by 2.3% and
%  the log-likelihood by 0.049, from fminbnd's relative tolerance on the
%  raw shift); the upper bound scaled by median(x - min(x)) instead of
%  range(x); both bounds evaluated exactly as candidates; a SmallSample
%  warning below n = 15.

    %% Constants and default settings

    gapFraction       = 0.001;
    minimumSampleSize = 5;
    smallSampleSize   = 15;   % Below this the min(x) spike often wins
    nGridPoints       = 200;  % Coarse scan resolution (in log gap)
    maxSeeds          = 5;    % Maximum number of local minima to refine

    % TolX is on t = log(min(x) - shift), which is dimensionless, so one
    % tolerance serves every unit and offset of x.
    optimizerOptions = optimset( ...
        'Display', 'off', ...
        'TolX',    1e-10, ...
        'TolFun',  1e-12, ...
        'MaxIter', 500);

    %% Validate and parse mode

    if nargin < 2 || isempty(mode)
        mode = '';
    end

    if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
        error('shiftlognormal_MLE:InvalidMode', ...
            'mode must be a character vector or string scalar.');
    end

    mode    = upper(string(mode));
    verbose = contains(mode, "VERBOSE") || contains(mode, "PLOT");

    %% Validate and clean the data

    validateattributes(x, {'numeric'}, ...
        {'real', 'vector', 'nonempty'}, mfilename, 'x', 1);

    x = x(:);

    finiteMask = isfinite(x);
    nRemoved   = sum(~finiteMask);
    x          = x(finiteMask);

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

    min_x     = min(x);
    max_x     = max(x);
    dataRange = max_x - min_x;

    if ~isfinite(dataRange) || dataRange <= 0
        error('shiftlognormal_MLE:ConstantData', ...
            ['The finite observations must have a positive range. ' ...
             'A shifted log-normal distribution cannot be fitted to ' ...
             'constant data.']);
    end

    if nSamples < smallSampleSize
        warning('shiftlognormal_MLE:SmallSample', ...
            ['Only %d observations. Below %d the unbounded spike of the ' ...
             'likelihood at shift -> min(x) often beats any interior ' ...
             'optimum (measured: 134/200 fits at n = 5, 9/200 at ' ...
             'n = 10, 0/200 at n = 15), so the shift is likely to sit ' ...
             'on its upper bound.'], nSamples, smallSampleSize);
    end

    %% Construct the data-derived shift bounds
    %
    % Work with offsets from the minimum, d = x - min(x) >= 0, and the gap
    % g = min(x) - shift > 0, so the shifted data are d + g. Forming them
    % this way (rather than x - shift) keeps full precision when x carries
    % a large offset, and searching over log(g) makes the fit equivariant
    % under x -> a*x + b.

    offsets    = x - min_x;
    gapScale   = median(offsets);       % typical distance above min(x)
    if ~(gapScale > 0)
        gapScale = dataRange;           % > half the data tied at min(x)
    end
    minimumGap  = gapFraction * gapScale;
    maximumGap  = dataRange;
    lowerShift  = min_x - maximumGap;
    upperShift  = min_x - minimumGap;
    shiftBounds = [lowerShift, upperShift];
    tBounds     = log([minimumGap, maximumGap]);   % t = log(gap)

    if verbose
        fprintf('Fitting a constrained shifted log-normal model...\n');
        fprintf('  n = %d,  data range = %.6g\n', ...
            nSamples, dataRange);
        fprintf('  Shift bounds: [%.8g, %.8g]\n', ...
            lowerShift, upperShift);
        fprintf('  Lower bound = min(x) - range(x)  [heuristic]\n');
        fprintf(['  Upper bound = min(x) - %.4g*median(x - min(x))  ' ...
                 '[gapFraction]\n'], ...
            gapFraction);
    end

    %% Coarse grid scan in t = log(gap) to find minimum basins
    %
    % shift = min(x) - exp(t). A grid uniform in t is dense near min(x),
    % where the profile changes fastest, and coarse far to the left, where
    % it flattens toward the Normal limit.

    objectiveFunction = @(t) neg_profile_loglike(offsets, exp(t));

    gridT   = linspace(tBounds(1), tBounds(2), nGridPoints);
    gridNLL = arrayfun(objectiveFunction, gridT);

    % Find interior local minima on the grid using strict inequality.
    isInteriorMin = [false, ...
        gridNLL(2:end-1) < gridNLL(1:end-2) & ...
        gridNLL(2:end-1) < gridNLL(3:end), ...
        false];

    % Always seed from the global grid minimum and add any local minima.
    [~, globalMinIdx] = min(gridNLL);
    seedIndices       = unique([globalMinIdx, find(isInteriorMin)]);

    % Sort seeds by grid NLL value and keep only the best maxSeeds.
    [~, sortOrd] = sort(gridNLL(seedIndices));
    seedIndices  = seedIndices(sortOrd(1:min(end, maxSeeds)));

    %% Refine each seed with fminbnd and select the best candidate
    %
    % Selection policy:
    %   1. If one or more refinements converge, return the converged
    %      candidate with the lowest NLL.
    %   2. If no refinement converges, return the finite nonconverged
    %      candidate with the lowest NLL and set Converged to false.
    %   3. The two bounds themselves are candidates too: fminbnd never
    %      evaluates the ends of its interval, so a solution on a bound
    %      would otherwise be returned slightly inside it, short of the
    %      optimum. A bound is taken only when its NLL is strictly lower
    %      than every refinement's; it then counts as converged, with
    %      ExitFlag 1 and an OptimizerOutput saying so.
    %
    % Consequently, Converged, ExitFlag, and OptimizerOutput always
    % describe the step that produced the returned solution.

    % Search radius spans two grid steps on each side of each seed.
    searchRadius = 2 * (tBounds(2) - tBounds(1)) / (nGridPoints - 1);

    % Track the best converged candidate separately from the best finite
    % candidate. The latter is used only if no refinement converges.
    bestConvergedT  = NaN;
    bestConvergedNLL    = Inf;
    bestConvergedFlag   = -1;
    bestConvergedOutput = struct( ...
        'message', 'No converged refinement completed.');

    bestFiniteT  = NaN;
    bestFiniteNLL    = Inf;
    bestFiniteFlag   = -1;
    bestFiniteOutput = struct( ...
        'message', 'No finite refinement completed.');

    for k = 1:numel(seedIndices)
        idx     = seedIndices(k);
        localLo = max(tBounds(1), gridT(idx) - searchRadius);
        localHi = min(tBounds(2), gridT(idx) + searchRadius);

        if localLo >= localHi
            continue;
        end

        [tK, nllK, flagK, outK] = fminbnd( ...
            objectiveFunction, localLo, localHi, optimizerOptions);

        candidateIsFinite = ...
            isfinite(tK) && isfinite(nllK);

        candidateConverged = ...
            candidateIsFinite && flagK > 0;

        % Track the lowest-NLL finite candidate regardless of convergence.
        % This candidate is used only if no refinement converges.
        if candidateIsFinite && nllK < bestFiniteNLL
            bestFiniteT  = tK;
            bestFiniteNLL    = nllK;
            bestFiniteFlag   = flagK;
            bestFiniteOutput = outK;
        end

        % Separately track the lowest-NLL converged candidate.
        if candidateConverged && nllK < bestConvergedNLL
            bestConvergedT  = tK;
            bestConvergedNLL    = nllK;
            bestConvergedFlag   = flagK;
            bestConvergedOutput = outK;
        end
    end

    % The bounds, evaluated exactly (policy 3).
    for tEnd = tBounds
        nllEnd = objectiveFunction(tEnd);
        if isfinite(nllEnd) && nllEnd < min(bestConvergedNLL, bestFiniteNLL)
            bestConvergedT      = tEnd;
            bestConvergedNLL    = nllEnd;
            bestConvergedFlag   = 1;
            bestConvergedOutput = struct('message', ...
                'Bound evaluated directly; it beat every refinement.');
        end
    end

    %% Select the returned solution

    if isfinite(bestConvergedT) && isfinite(bestConvergedNLL)
        % At least one refinement converged. Return the converged candidate
        % with the lowest NLL, even if a nonconverged candidate reported a
        % slightly lower objective value.
        bestT           = bestConvergedT;
        minNLL          = bestConvergedNLL;
        exitFlag        = bestConvergedFlag;
        optimizerOutput = bestConvergedOutput;
        converged       = true;

    elseif isfinite(bestFiniteT) && isfinite(bestFiniteNLL)
        % No refinement converged, but at least one finite candidate was
        % obtained. Return the finite candidate with the lowest NLL.
        bestT           = bestFiniteT;
        minNLL          = bestFiniteNLL;
        exitFlag        = bestFiniteFlag;
        optimizerOutput = bestFiniteOutput;
        converged       = false;

    else
        % No finite candidate was obtained.
        error('shiftlognormal_MLE:OptimizationFailed', ...
            ['All %d refinement attempts returned nonfinite results. ' ...
             'The optimization failed completely.'], ...
            numel(seedIndices));
    end

    bestGap   = exp(bestT);
    bestShift = min_x - bestGap;

    if verbose
        fprintf(['  Grid seeds refined: %d  |  Best shift: %.8g  |  ' ...
                 'NLL: %.6g  |  Converged: %d\n'], ...
            numel(seedIndices), bestShift, minNLL, converged);
    end

    %% Warn if the returned solution did not converge

    if ~converged
        warning('shiftlognormal_MLE:PartialConvergence', ...
            ['No refinement pass converged. The finite candidate with ' ...
             'the lowest NLL is being returned (exit flag %d, ' ...
             'message: %s). ModelResults.Converged is false. Inspect ' ...
             'ModelResults.Diagnostics and consider adjusting the shift ' ...
             'bounds or optimizer options.'], ...
            exitFlag, get_optimizer_message(optimizerOutput));
    end

    %% Recover mu and sigma at the fitted shift

    shiftedData = offsets + bestGap;   % = x - bestShift, without cancellation

    if any(shiftedData <= 0) || any(~isfinite(shiftedData))
        error('shiftlognormal_MLE:InvalidShiftedData', ...
            ['The fitted shift produced nonpositive or nonfinite ' ...
             'shifted observations.']);
    end

    logData = log(shiftedData);
    bestMu  = mean(logData);

    % std(..., 1) divides by n, giving the MLE of sigma.
    bestSigma = std(logData, 1);

    if ~isfinite(bestMu) || ~isfinite(bestSigma) || bestSigma <= 0
        error('shiftlognormal_MLE:InvalidParameterEstimate', ...
            'The fitted mu or sigma estimate is invalid.');
    end

    %% Calculate likelihood and information criteria

    logLikelihood = -minNLL;
    nParameters   = 3;

    % AIC = 2k - 2*log(L)
    AIC = 2*nParameters - 2*logLikelihood;

    % BIC = k*log(n) - 2*log(L)
    BIC = nParameters*log(nSamples) - 2*logLikelihood;

    %% Diagnose proximity to shift bounds

    % Judged in t = log(gap), which is dimensionless, so the same tolerance
    % applies at every unit and offset of x. BoundaryTolerance is in t.

    boundaryTolerance = 1e-6 * (tBounds(2) - tBounds(1));

    distanceToLowerBound = bestShift - lowerShift;
    distanceToUpperBound = upperShift - bestShift;

    nearLowerBound = ...
        tBounds(2) - bestT <= boundaryTolerance;

    nearUpperBound = ...
        bestT - tBounds(1) <= boundaryTolerance;

    nearAnyBoundary = nearLowerBound || nearUpperBound;
    boundaryMessage = '';

    if nearLowerBound
        boundaryMessage = sprintf( ...
            ['The estimated shift (%.8g) is at or near the lower bound ' ...
             '(%.8g = min(x) - range(x)). The true shift may lie further ' ...
             'left. The result is clipped and depends on how this bound ' ...
             'was defined. Consider widening the lower bound.'], ...
            bestShift, lowerShift);

        warning('shiftlognormal_MLE:ShiftNearLowerBound', ...
            '%s', boundaryMessage);

    elseif nearUpperBound
        boundaryMessage = sprintf( ...
            ['The estimated shift (%.8g) is at or near the upper bound ' ...
             '(%.8g = min(x) - %.4g*median(x - min(x))). This is the ' ...
             'unbounded spike of the likelihood at min(x), not an ' ...
             'interior optimum: the result depends on gapFraction. ' ...
             'Common below n = %d; otherwise check the data.'], ...
            bestShift, upperShift, gapFraction, smallSampleSize);

        warning('shiftlognormal_MLE:ShiftNearUpperBound', ...
            '%s', boundaryMessage);
    end

    %% Calculate fitted PDF and corresponding coordinates

    [fit_x, sortIdx] = sort(x);

    % Evaluated on the shifted values formed without cancellation.
    fitLogPDF = shifted_lognormal_logpdf( ...
        shiftedData(sortIdx), bestMu, bestSigma);

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

    ModelResults.Mu    = bestMu;
    ModelResults.Sigma = bestSigma;
    ModelResults.Shift = bestShift;

    ModelResults.LogLike = logLikelihood;
    ModelResults.AIC     = AIC;
    ModelResults.BIC     = BIC;

    ModelResults.FitX = fit_x;
    ModelResults.Fit  = fitPDF;

    ModelResults.N           = nSamples;
    ModelResults.NParameters = nParameters;

    % These three fields all describe the refinement that produced the
    % returned parameter estimates.
    ModelResults.Converged       = converged;
    ModelResults.ExitFlag        = exitFlag;
    ModelResults.OptimizerOutput = optimizerOutput;

    %% Store shift and boundary diagnostics

    ModelResults.Diagnostics = struct();

    ModelResults.Diagnostics.MinimumObservation = min_x;
    ModelResults.Diagnostics.MaximumObservation = max_x;
    ModelResults.Diagnostics.DataRange          = dataRange;

    ModelResults.Diagnostics.GapFraction = gapFraction;
    ModelResults.Diagnostics.MinimumGap  = minimumGap;

    ModelResults.Diagnostics.ShiftBounds = shiftBounds;
    ModelResults.Diagnostics.LowerShift  = lowerShift;
    ModelResults.Diagnostics.UpperShift  = upperShift;

    ModelResults.Diagnostics.LowerBoundNote = ...
        ['lowerShift = min(x) - range(x).  Heuristic: shift clipped ' ...
         'if true value lies further left.'];

    ModelResults.Diagnostics.UpperBoundNote = sprintf( ...
        ['upperShift = min(x) - %.4g*median(x - min(x)).  Controlled ' ...
         'by gapFraction (range(x) replaces the median when more than ' ...
         'half the data tie at min(x)).'], ...
        gapFraction);

    ModelResults.Diagnostics.GapScale = gapScale;

    ModelResults.Diagnostics.ShiftToMinimumGap = ...
        bestGap;

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

    ModelResults.Diagnostics.GridPoints = ...
        nGridPoints;

    ModelResults.Diagnostics.GridSeedsRefined = ...
        numel(seedIndices);

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

    function nll = neg_profile_loglike(offsetsIn, gap)
        % NEG_PROFILE_LOGLIKE Profile negative log-likelihood over the gap.
        %
        % With x - shift = (x - min(x)) + gap, mu and sigma are
        % concentrated out analytically:
        %
        %   NLL = n*log(sigma_hat) + sum(log(x - shift)) ...
        %       + n/2*(1 + log(2*pi))
        %
        % A large finite penalty is returned for degenerate inputs so that
        % fminbnd receives a continuous objective without hard Inf jumps.

        shiftedValues = offsetsIn + gap;

        if any(shiftedValues <= 0) || any(~isfinite(shiftedValues))
            nll = 1e15;
            return;
        end

        logShifted = log(shiftedValues);
        sigmaMLE   = std(logShifted, 1);

        if ~isfinite(sigmaMLE) || sigmaMLE <= 0
            nll = 1e15;
            return;
        end

        n = numel(offsetsIn);

        nll = n*log(sigmaMLE) ...
            + sum(logShifted) ...
            + 0.5*n*(1 + log(2*pi));

        if ~isfinite(nll)
            nll = 1e15;
        end
    end

    function logPDF = shifted_lognormal_logpdf(shiftedZ, mu, sigma)
        % SHIFTED_LOGNORMAL_LOGPDF Log-PDF of the shifted log-normal,
        % taking the already-shifted values z - shift.
        %
        % Points at or below the shift receive log-density -Inf.

        logPDF   = -Inf(size(shiftedZ));
        valid    = shiftedZ > 0 & isfinite(shiftedZ);

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
        % GET_OPTIMIZER_MESSAGE Retrieve the fminbnd termination message.

        if isstruct(output) && ...
                isfield(output, 'message') && ...
                ~isempty(output.message)

            message = output.message;
        else
            message = 'No optimizer message was returned.';
        end
    end
end
