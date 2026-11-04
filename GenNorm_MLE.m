function ModelResults = GenNorm_MLE(X, models_to_test, mode)
% GENNORM_MLE Fits a Generalized Normal (Gaussian) Distribution mixture
% model.
%
%   ModelResults = GENNORM_MLE(X, MODELS_TO_TEST, MODE)
%
%   The probability density function (PDF) for a single term is:
%       f(x) = alpha / (2*sigma*Gamma(1/alpha)) * exp( -
%       (|x-mu|/sigma)^alpha )
%
%   INPUTS:
%       X              - Numeric data vector. MODELS_TO_TEST - Vector of
%       component counts (e.g., [1, 2, 3]). MODE           - Flags: 'Fit',
%       'Plot', 'Select'.
%
%   OUTPUT:
%       ModelResults   - Structure array with Params [mu, sigma, alpha,
%       weight].

%% 1. Input Sanitization
if nargin < 2 || isempty(models_to_test)
    models_to_test = 1;
    fprintf('Default: Fitting 1 term.\n');
end
if nargin < 3
    mode = '';
end

X = X(:);
X = X(~isnan(X) & ~isinf(X));
if isempty(X)
    error('Input data must be numeric and non-empty.');
end

models_to_test = sort(unique(models_to_test));

% Initialize Output Structure
ModelResults = repmat(struct('Model',[],'Params',[],'LogLike',[],...
    'AIC',[], 'BIC',[], 'Fit',[]), ...
    numel(models_to_test), 1);

%% 2. Main Optimization Loop
for m = 1:length(models_to_test)
    kComp = models_to_test(m);

    fprintf('Fitting %d-term Generalized Normal model...\n', kComp);

    % --- Run the Engine ---
    stats = FitGenNorm_Engine(X, kComp);

    % --- Store Results ---
    ModelResults(m).Model    = sprintf('%d-term Generalized Normal (mu, sigma, alpha, weight)', kComp);
    ModelResults(m).Params   = stats.Params; % [mu, sigma, alpha, weight]
    ModelResults(m).LogLike  = stats.LogLike;

    % --- Info Criteria ---
    nParams = 4 * kComp - 1; % (mu, sig, alpha) per comp + (k-1) weights
    nObs = numel(X);
    LL = stats.LogLike;

    ModelResults(m).AIC = 2*nParams - 2*LL + (2*nParams*(nParams+1))/(nObs-nParams-1);
    ModelResults(m).BIC = nParams*log(nObs) - 2*LL;
end

%% 3. Model Selection (Likelihood Ratio Test)
if contains(upper(mode),'SELECT') && length(models_to_test) > 1
    % Greedy forward selection
    best_idx = 1;
    alpha_crit = 0.01;

    for i = 1:length(models_to_test)-1
        % Null: Simpler Model (i), Alt: Complex Model (i+1)
        LL_null = ModelResults(i).LogLike;
        LL_alt  = ModelResults(i+1).LogLike;

        df = 4 * (models_to_test(i+1) - models_to_test(i));
        LR_stat = max(0, 2 * (LL_alt - LL_null));
        pval = 1 - chi2cdf(LR_stat, df);

        if pval < alpha_crit
            best_idx = i + 1; % Significant improvement, keep expanding
        else
            break; % No improvement, stop here
        end
    end
    ModelResults = ModelResults(best_idx);
    fprintf('>> Auto-Selected: %s (p=%.4g)\n', ModelResults.Model, pval);
end

%% 4. Plotting
if contains(upper(mode),'PLOT')
    plot_results(X, ModelResults);
end

%% 5. Return Fit Curve
if contains(upper(mode),'FIT')
    sorted_x = sort(X);
    for i = 1:numel(ModelResults)
        ModelResults(i).Fit = EvaluateGenNorm(sorted_x, ModelResults(i).Params);
    end
end
end

%% ============================================================
%  CORE OPTIMIZATION ENGINE
% ============================================================
function stats = FitGenNorm_Engine(X, kComp)

% --- 1. Define Bounds & Constants --- Sigma Floor: Based on IQR to avoid
% outliers causing singularities
range_iqr = iqr(X);
if range_iqr == 0, range_iqr = std(X); end
sig_min = max(1e-6, 1e-3 * range_iqr);

% Logit Bounds (Weights)
max_logit = 20;

% Parameter Packing: P = [Mu(1:k), LogSigma(1:k), LogAlpha(1:k),
% Logits(1:k-1)]

% Bounds for FMINCON Mu: Unconstrained (effectively) LogSigma: >
% log(sig_min) LogAlpha: log(0.5) to log(10) (Standard GenNorm range)

lb_mu    = -inf(1, kComp);
lb_sig   = log(sig_min) * ones(1, kComp);
lb_alpha = log(0.5) * ones(1, kComp);     % Alpha >= 0.5 (very peaked)
lb_w     = -max_logit * ones(1, kComp-1);

ub_mu    = inf(1, kComp);
ub_sig   = inf(1, kComp);
ub_alpha = log(10) * ones(1, kComp);      % Alpha <= 10 (very boxy)
ub_w     = max_logit * ones(1, kComp-1);

LB = [lb_mu, lb_sig, lb_alpha, lb_w];
UB = [ub_mu, ub_sig, ub_alpha, ub_w];

% --- 2. Generate Robust Seeds ---
seeds = GenerateSeeds(X, kComp, sig_min, max_logit);

% --- 3. Optimization Settings ---
opt = optimoptions('fmincon', ...
    'Display','off', ...
    'Algorithm','interior-point', ... % Robust to bounds
    'MaxIterations', 2000, ...
    'OptimalityTolerance', 1e-5);

bestLL = -Inf;
bestP  = [];

% --- 4. Multi-Start Loop ---
nStarts = size(seeds, 1);
warning('off', 'MATLAB:nearlySingularMatrix');

for i = 1:nStarts
    p0 = seeds(i,:);

    % Objective Function
    obj = @(p) GenNorm_NLL(X, p, kComp);

    try
        [popt, nll] = fmincon(obj, p0, [], [], [], [], LB, UB, [], opt);

        if -nll > bestLL
            bestLL = -nll;
            bestP  = popt;
        end
    catch
        % Skip failed starts
    end
end
warning('on', 'MATLAB:nearlySingularMatrix');

% --- 5. Reconstruct Parameters ---
if isinf(bestLL)
    stats.Params  = NaN(kComp, 4);
    stats.LogLike = -Inf;
else
    stats.Params  = ReconstructParams(bestP, kComp);
    stats.LogLike = bestLL;
end
end

% ============================================================
%  LIKELIHOOD FUNCTION (Vectorized & Stable)
% ============================================================
function nll = GenNorm_NLL(x, p, kComp)
% Unpack Parameters
mu       = p(1:kComp);
sigma    = exp(p(kComp+1:2*kComp));
alpha    = exp(p(2*kComp+1:3*kComp));

if kComp > 1
    v = [p(3*kComp+1:end), 0];
    v = v - max(v); % Shift for stability
    w = exp(v) / sum(exp(v));
else
    w = 1;
end

% --- GenNorm Log-PDF Calculation --- log f(x) = log(alpha) -
% log(2*sigma*Gamma(1/alpha)) - (|x-mu|/sigma)^alpha
%          = log(alpha) - log(2) - log(sigma) - gammaln(1/alpha) - ...

n = numel(x);
log_w = log(w + 1e-100);

% Precompute constant terms for each component (1 x k)
log_const = log(alpha) - log(2) - log(sigma) - gammaln(1./alpha);

% Compute variable term (n x k) Matrix expansion: (x - mu)
log_pdf = zeros(n, kComp);

for k = 1:kComp
    z = abs(x - mu(k)) ./ sigma(k);
    log_pdf(:, k) = log_w(k) + log_const(k) - (z .^ alpha(k));
end

% Log-Sum-Exp Trick
max_L = max(log_pdf, [], 2);
sum_exp = sum(exp(log_pdf - max_L), 2);
log_L = max_L + log(sum_exp);

nll = -sum(log_L);
if isnan(nll), nll = Inf; end
end

% ============================================================
%  SEED GENERATION (Shape Scanning)
% ============================================================
function seeds = GenerateSeeds(x, kComp, sig_min, max_logit)
% Generates guesses covering different "Shape Regimes": 1. Laplacian
% (alpha=1) 2. Gaussian (alpha=2) 3. Platykurtic (alpha=4)

seeds = [];

mu_global = mean(x);
std_global = std(x);

% --- Strategy 1: K-Means (Standard) --- We run K-means to get locations,
% then permute Alpha
warning('off','stats:kmeans:FailedToConverge');
try
    [idx, C] = kmeans(x, kComp, 'Replicates', 3, 'MaxIter', 100, 'Display','off');
    [mu_k, ord] = sort(C'); %#ok<*TRSRT>

    sig_k = zeros(1, kComp);
    w_k   = zeros(1, kComp);

    for j = 1:kComp
        cluster_data = x(idx == ord(j));
        if numel(cluster_data) > 1
            sig_k(j) = std(cluster_data);
        else
            sig_k(j) = std_global/kComp;
        end
        w_k(j) = numel(cluster_data)/numel(x);
    end
    sig_k = max(sig_k, sig_min);
    w_k = max(w_k, 0.01); w_k = w_k/sum(w_k);

    % Convert weights to logits
    if kComp > 1
        v_k = log(w_k(1:end-1)) - log(w_k(end));
        v_k = max(min(v_k, max_logit), -max_logit);
    else
        v_k = [];
    end

    % Create 3 Variants of this seed (Laplace, Gauss, Box)
    alphas = [1.0, 2.0, 5.0];
    for a_val = alphas
        row = [mu_k, log(sig_k), log(repmat(a_val,1,kComp)), v_k];
        seeds = [seeds; row]; %#ok<AGROW>
    end

catch
    % Fallback if Kmeans fails
end
warning('on','stats:kmeans:FailedToConverge');

% --- Strategy 2: Quantiles (Spread) ---
if kComp > 1
    mu_q = quantile(x, linspace(0.1, 0.9, kComp));
    sig_q = repmat(std_global/kComp, 1, kComp);
    if kComp > 1, v_equal = zeros(1, kComp-1); else, v_equal = []; end

    seeds = [seeds; mu_q, log(sig_q), log(2*ones(1,kComp)), v_equal];
end

% --- Strategy 3: Global (One Big Hump) --- Good for "Peaked" data centered
% at 0
mu_g = repmat(mu_global, 1, kComp) + linspace(-1,1,kComp)*std_global*0.1;
sig_g = repmat(std_global, 1, kComp);
if kComp > 1, v_equal = zeros(1, kComp-1); else, v_equal = []; end

seeds = [seeds; sort(mu_g), log(sig_g), log(1.5*ones(1,kComp)), v_equal];

% Safety
if isempty(seeds)
    seeds = [linspace(min(x),max(x),kComp), log(repmat(std(x),1,kComp)), log(2*ones(1,kComp)), zeros(1,max(0,kComp-1))];
end
end

% ============================================================
%  UTILITIES
% ============================================================
function Params = ReconstructParams(p, kComp)
mu    = p(1:kComp);
sigma = exp(p(kComp+1:2*kComp));
alpha = exp(p(2*kComp+1:3*kComp));

if kComp > 1
    v = [p(3*kComp+1:end), 0];
    v = v - max(v);
    w = exp(v) / sum(exp(v));
else
    w = 1;
end

% Return Nx4 Matrix
Table = [mu(:), sigma(:), alpha(:), w(:)];
Params = sortrows(Table, 1);
end

function y = EvaluateGenNorm(x, Params)
y = zeros(size(x));
for k = 1:size(Params,1)
    mu  = Params(k,1);
    sig = Params(k,2);
    alp = Params(k,3);
    w   = Params(k,4);

    % PDF Calculation
    coef = alp / (2*sig*gamma(1/alp));
    term = exp( - (abs(x-mu)./sig).^alp );
    y = y + w * coef * term;
end
end

function plot_results(x, res)
figure('Color','w'); hold on;
histogram(x, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.2 0.2 0.2]);

xx = linspace(min(x), max(x), 500)';
colors = lines(numel(res));

for i = 1:numel(res)
    yy = EvaluateGenNorm(xx, res(i).Params);
    plot(xx, yy, 'LineWidth', 2, 'Color', colors(i,:), ...
        'DisplayName', res(i).Model);
end
legend('show'); grid on; xlabel('Data'); ylabel('PDF');
title('Generalized Normal Fit');
end