function ModelResults = shiftweibull_MLE(x, mode)
% SHIFTWEIBULL_MLE Maximum Likelihood for 3-Parameter Weibull.
%
%   ModelResults = shiftweibull_MLE(x, mode)
%
%   PDF: f(x) = (b/a) * ((x-c)/a)^(b-1) * exp( -((x-c)/a)^b )
%
%   Improvements in this version:
%   - Simultaneous Optimization: Eliminates the slow nested wblfit loop.
%   - Unconstrained Mapping: Enforces b >= 1.001 to avoid the degenerate 
%     unbounded likelihood trap common in 3-parameter Weibulls.
%   - Bootstrap SE & CI: Computes standard errors robustly, avoiding the 
%     Hessian singularity issues when shape (b) <= 2.

%% 1. Input Handling
if nargin < 2, mode = ''; end
x = x(:);
x = x(~isnan(x) & ~isinf(x));
n = length(x);

if n < 4
    error('At least 4 valid data points are required.');
end

ModelResults = struct('Model', 'Shifted-Weibull (Scale, Shape, Shift)', ...
                      'Params', [], 'LogLike', [], 'AIC', [], 'BIC', [], ...
                      'SE', [], 'CILower', [], 'CIUpper', [], ...
                      'CIMethod', 'Nonparametric Bootstrap (n=300)', ...
                      'Fit', [], 'FitX', []);

fprintf('Fitting Shifted Weibull model...\n');

%% 2. Robust Initialization & Simultaneous Optimization
min_x = min(x);

% Initial Guess: Shift slightly below min(x)
c_init = min_x - 0.05 * std(x);
shifted_data = x - c_init;

% Quick heuristic for Scale (a) and Shape (b)
try
    par_init = wblfit(shifted_data);
    a_init = par_init(1);
    b_init = max(par_init(2), 1.1); % Ensure safely > 1
catch
    a_init = mean(shifted_data);
    b_init = 2.0;
end

% Map to unconstrained space (q1, q2, q3)
q0 = [log(a_init), log(b_init - 1.001), log(min_x - c_init)];

% Optimize all 3 parameters simultaneously (Analytical Vectorized NLL)
opt = optimset('Display', 'off', 'MaxFunEvals', 3000, 'TolX', 1e-6);
obj = @(q) nll_shiftweibull(q, x, min_x);

[q_best, nll_min] = fminsearch(obj, q0, opt);

% Recover optimal physical parameters
a_hat = exp(q_best(1));
b_hat = 1.001 + exp(q_best(2));
c_hat = min_x - exp(q_best(3));

ModelResults.Params  = [a_hat, b_hat, c_hat];
ModelResults.LogLike = -nll_min;

% Info Criteria
k_params = 3; 
ModelResults.AIC = 2*k_params + 2*nll_min + (2*k_params*(k_params+1))/(n-k_params-1);
ModelResults.BIC = k_params*log(n) + 2*nll_min;

%% 3. Robust Standard Errors & CI via Nonparametric Bootstrap
fprintf('Calculating Robust Bootstrap CIs (may take a few seconds)...\n');
nBoot = 300; 
boot_params = zeros(nBoot, 3);

% Turn off warnings for bootstrap loop
warning('off', 'all'); 
for i = 1:nBoot
    % Resample with replacement
    xb = x(randi(n, n, 1));
    min_xb = min(xb);
    
    % Adjust starting q3 if resampled min(x) changed drastically
    q_boot_init = q_best;
    if min_xb <= c_hat
        q_boot_init(3) = log(min_xb - (min_xb - 0.01*std(xb)));
    end
    
    obj_b = @(q) nll_shiftweibull(q, xb, min_xb);
    q_b = fminsearch(obj_b, q_boot_init, opt);
    
    boot_params(i, 1) = exp(q_b(1));
    boot_params(i, 2) = 1.001 + exp(q_b(2));
    boot_params(i, 3) = min_xb - exp(q_b(3));
end
warning('on', 'all');

% Calculate SE and 95% Confidence Intervals
ModelResults.SE      = std(boot_params, 0, 1);
ModelResults.CILower = prctile(boot_params, 2.5, 1);
ModelResults.CIUpper = prctile(boot_params, 97.5, 1);

%% 4. Plotting / Fitting
if contains(upper(mode), 'FIT') || contains(upper(mode), 'PLOT')
    xs = sort(x);
    
    % Analytic PDF computation (bypasses dependency on wblpdf)
    z = (xs - c_hat) ./ a_hat;
    y_fit = (b_hat / a_hat) .* (z.^(b_hat - 1)) .* exp(-(z.^b_hat));
    
    ModelResults.FitX = xs;
    ModelResults.Fit  = y_fit;
end

if contains(upper(mode), 'PLOT')
    figure('Color', 'w');
    histogram(x, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
              'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.5);
    hold on;
    plot(xs, y_fit, 'r-', 'LineWidth', 2);
    xlabel('Data Values'); ylabel('Probability Density');
    title(sprintf('Shifted Weibull Fit\nScale(a)=%.2f, Shape(b)=%.2f, Shift(c)=%.2f', ...
          a_hat, b_hat, c_hat));
    legend('Empirical PDF', 'MLE Fit', 'Location', 'best');
    grid on; hold off;
end

end

%% ============================================================
%  INTERNAL FULLY VECTORIZED OBJECTIVE FUNCTION
% ============================================================
function nll = nll_shiftweibull(q, x, min_x)
    % 1. Inverse Unconstrained Mapping
    a = exp(q(1));
    b = 1.001 + exp(q(2)); % Strictly bounds b >= 1.001
    c = min_x - exp(q(3)); % Strictly bounds c < min(x)
    
    % 2. Standardize data
    z = (x - c) ./ a;
    
    % 3. Vectorized Log-Likelihood Calculation
    % log(f) = log(b) - log(a) + (b-1)*log(z) - z^b
    log_pdf = log(b) - log(a) + (b - 1) .* log(z) - (z.^b);
    
    nll = -sum(log_pdf);
    
    % Safety catch for numerical precision limits
    if ~isfinite(nll)
        nll = 1e99; 
    end
end