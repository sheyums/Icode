function ModelResults = multigaussian_MLE(x, testmodels, mode)
% MULTIGAUSSIAN_MLE  Robust Maximum Likelihood Estimation of Gaussian
% Mixtures
%
%   ModelResults = multigaussian_MLE(x, testmodels, mode)
%
%   Fits Gaussian mixture models (GMMs) with varying numbers of components
%   using direct maximum likelihood optimization with multiple starts.
%
%   Features:
%     - Numerically stable log-sum-exp likelihood - Softmax
%     parameterization of weights (no constraints) - Multi-start
%     deterministic seeding (K-means + quantiles) - BIC-based model
%     selection (recommended for mixtures) - Optional plotting and curve
%     output
%
%   INPUTS:
%     x          : data vector testmodels : vector of component counts
%     (e.g. [1 2 3]) mode       : flags string:
%                   'Select' → choose best model by BIC 'Plot'   → plot
%                   histogram + fits 'Fit'    → return fitted pdf values
%
%   OUTPUT:
%     ModelResults : structure array (or scalar if Select enabled)


%% -------------------- Input handling --------------------
if nargin < 2 || isempty(testmodels)
    testmodels = 1;
end
if nargin < 3
    mode = '';
end

x = x(:);
x = x(~isnan(x));
x = x(~isinf(x)); % Safety: remove Infs
n = numel(x);

if isempty(x)
    error('Input data vector is empty after removing NaNs/Infs.');
end

testmodels = sort(unique(testmodels(:)')); % Ensure unique and sorted
ModelResults = repmat(struct( ...
    'Model', [], ...
    'Params', [], ...
    'LogLike', [], ...
    'AIC', [], ...
    'BIC', [], ...
    'Fit', []), numel(testmodels), 1);

%% -------------------- Main loop --------------------
for m = 1:numel(testmodels)
    kComp = testmodels(m);

    if kComp == 1
        fprintf('Fitting model with %d Gaussian term...\n', kComp);
        % --- Analytic Shortcut for k=1 (Instant & Exact) ---
        mu_hat  = mean(x);
        sig_hat = std(x);
        %logL    = sum(log(normpdf(x, mu_hat, sig_hat)));
        logL    = -normlike([mu_hat, sig_hat],x);


        ModelResults(m).Model    = '1-term Gaussian (mu, sigma)';
        ModelResults(m).Params   = [mu_hat, sig_hat, 1];
        ModelResults(m).LogLike  = logL;
    else
        % --- Numerical Optimization for k > 1 ---
        fprintf('Fitting model with %d Gaussian term(s)...\n', kComp);
        stats = FitGMM_Engine(x, kComp);

        ModelResults(m).Model    = [num2str(kComp),'-term Gaussian (mu, sigma, weight)'];
        ModelResults(m).Params   = stats.Params;
        ModelResults(m).LogLike  = stats.LogLike;
    end

    % --- Information Criteria ---
    k = 3*kComp - 1; % Free parameters (mu + sigma + weights-1)
    LL = ModelResults(m).LogLike;

    ModelResults(m).AIC = 2*k - 2*LL + (2*k*(k+1))/(n-k-1);
    ModelResults(m).BIC = k*log(n) - 2*LL;
end

%% -------------------- Model selection (BIC) --------------------
if contains(upper(mode), 'SELECT') && numel(ModelResults) > 1
    [~, idx] = min([ModelResults.BIC]);
    bestModel = ModelResults(idx);
    fprintf('>> Selected %s based on BIC.\n', bestModel.Model);
    ModelResults = bestModel;
end

%% -------------------- Plotting --------------------
if contains(upper(mode), 'PLOT')
    % Create a specialized visualization
    figure('Color','w'); hold on;

    % 1. Data Histogram
    histogram(x, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
        'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5);

    % 2. Smooth evaluation grid
    margin = (max(x) - min(x)) * 0.2;
    xs = linspace(min(x)-margin, max(x)+margin, 1000)';

    colors = lines(numel(ModelResults));

    % 3. Plot Models
    for i = 1:numel(ModelResults)
        y = EvaluateGMM(xs, ModelResults(i).Params);
        plot(xs, y, 'LineWidth', 2, 'Color', colors(i,:), ...
            'DisplayName', sprintf('k=%d (BIC=%.1f)', size(ModelResults(i).Params,1), ModelResults(i).BIC));
    end

    xlabel('Data Value'); ylabel('Probability Density');
    title('Gaussian Mixture Model Fit');
    legend('show', 'Location', 'best');
    grid on; box on; hold off;
end

%% -------------------- Optional fitted curve --------------------
if contains(upper(mode), 'FIT')
    % Compute PDF at specific input points (sorted for clean plotting by
    % user)
    xs = sort(x);
    for i = 1:numel(ModelResults)
        ModelResults(i).Fit = EvaluateGMM(xs, ModelResults(i).Params);
    end
end

end

%% =====================================================================
%                           INTERNAL FUNCTIONS
% =====================================================================

function stats = FitGMM_Engine(x, kComp)
% 1. Robust Sigma Floor (Based on IQR to ignore outliers) Prevents
% "infinite spike" singularities
data_iqr = iqr(x);
if data_iqr == 0, data_iqr = std(x); end
sig_min = max(1e-6, 1e-3 * data_iqr);

% 2. Optimization Bounds Logits constrained to +/- 20 (approx weight range
% 2e-9 to 1)
max_logit = 20;

% Parameter Packing: [Mu(1:k), LogSigma(1:k), Logits(1:k-1)] Note: LogSigma
% allows unconstrained optimization of scale

% LB: [Mu=-inf, LogSigma=log(min), Logits=-20]
LB = [-inf(1,kComp), log(sig_min)*ones(1,kComp), -max_logit*ones(1,kComp-1)];

% UB: [Mu=+inf, LogSigma=+inf,     Logits=+20]
UB = [ inf(1,kComp), inf(1,kComp),                max_logit*ones(1,kComp-1)];

% 3. Optimizer Settings 'interior-point' handles bounds robustly
opt = optimoptions('fmincon', ...
    'Display','off', ...
    'Algorithm','interior-point', ...
    'MaxIterations', 1000, ...
    'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-5);

% 4. Generate Robust Seeds
seeds = GenerateSeeds(x, kComp, sig_min, max_logit);
nStarts = size(seeds, 1);

all_LL = -inf(nStarts, 1);
all_P  = nan(nStarts, 3*kComp - 1);

% 5. Multi-start Optimization
warning('off', 'MATLAB:nearlySingularMatrix');

for i = 1:nStarts
    p0 = seeds(i,:);

    % Define Objective
    obj = @(p) GMM_NegLogLike(x, p, kComp);

    try
        [popt, nll] = fmincon(obj, p0, [], [], [], [], LB, UB, [], opt);
        all_LL(i)  = -nll;
        all_P(i,:) = popt;
    catch
        % If optimization crashes (rare with bounds), skip
    end
end

warning('on', 'MATLAB:nearlySingularMatrix');

% 6. Select Best Result
[bestLL, idx] = max(all_LL);

if ~isfinite(bestLL)
    % Fallback for total failure (return NaNs)
    stats.Params  = NaN(kComp, 3);
    stats.LogLike = -Inf;
else
    stats.Params  = ReconstructParams(all_P(idx,:), kComp);
    stats.LogLike = bestLL;
end
end

%% -------------------- Negative Log-Likelihood --------------------
function nll = GMM_NegLogLike(x, p, kComp)
% Extract Parameters
mu    = p(1:kComp);
sigma = exp(p(kComp+1:2*kComp)); % Enforces sigma > 0

if kComp > 1
    % Softmax weights: exp(v) / sum(exp(v)) We append 0 for the k-th
    % component (reference)
    v = [p(2*kComp+1:end), 0];

    % Numerical stability: subtract max(v) before exp
    v_shift = v - max(v);
    w = exp(v_shift) / sum(exp(v_shift));
else
    w = 1;
end

% Log-Likelihood Calculation
n = numel(x);
log_w = log(w + 1e-300); % Safety floor
log_sig = log(sigma);
c = -0.5 * log(2*pi);

% Compute log-probability for every point x and component k Matrix way:
% size(n, k)
log_pdf_k = zeros(n, kComp);

for k = 1:kComp
    z = (x - mu(k)) ./ sigma(k);
    log_pdf_k(:,k) = log_w(k) + c - log_sig(k) - 0.5 * z.^2;
end

% Log-Sum-Exp Trick for stability: log(sum(exp(vals))) log(sum(w_k * N(x)))
% = M + log(sum(exp(log_pdf_k - M)))
M = max(log_pdf_k, [], 2);
log_lik_per_point = M + log(sum(exp(log_pdf_k - M), 2));

nll = -sum(log_lik_per_point);

if isnan(nll), nll = Inf; end
end

%% -------------------- Robust Seed Generation --------------------
function seeds = GenerateSeeds(x, kComp, sig_min, max_logit)
nParams = 3*kComp - 1;
seeds = [];

% Global stats
mu_global = mean(x);
std_global = std(x);

% --- Strategy 1: K-means (Primary) --- Suppress convergence warnings; we
% only need rough starting points.
warning('off', 'stats:kmeans:FailedToConverge');

for r = 1:4
    try
        % Increased MaxIter to 500 to reduce warnings Added EmptyAction ->
        % singleton to prevent crashes on rare data
        [idx, C] = kmeans(x, kComp, 'Replicates', 1, ...
            'MaxIter', 500, ...
            'Display', 'off', ...
            'EmptyAction', 'singleton');

        [mu0, order] = sort(C'); %#ok<*TRSRT>

        sig0 = zeros(1, kComp);
        w0   = zeros(1, kComp);

        for j = 1:kComp
            cluster_data = x(idx == order(j));
            if numel(cluster_data) > 1
                sig0(j) = std(cluster_data);
            else
                sig0(j) = std_global / kComp;
            end
            w0(j) = numel(cluster_data) / numel(x);
        end

        % Enforce floors
        sig0 = max(sig0, sig_min);
        w0   = max(w0, 0.01);
        w0   = w0 / sum(w0);

        % Convert to params
        if kComp > 1
            v0 = log(w0(1:end-1)) - log(w0(end));
            v0 = max(min(v0, max_logit), -max_logit); % Clamp
        else
            v0 = [];
        end

        row_seed = [mu0(:)', log(sig0(:)'), v0(:)'];

        if length(row_seed) == nParams
            seeds = [seeds; row_seed]; %#ok<*AGROW>
        end
    catch
        % If K-means fails entirely, just skip this seed
    end
end

% Restore warnings
warning('on', 'stats:kmeans:FailedToConverge');

% --- Strategy 2: Uniform Spacing (Grid) ---
if kComp == 1
    mu_grid = mu_global;
else
    mu_grid = linspace(min(x), max(x), kComp + 2);
    mu_grid = mu_grid(2:end-1);
end

sig_grid = repmat((max(x)-min(x))/(2*kComp), 1, kComp);

if kComp > 1
    v_equal = zeros(1, kComp-1);
else
    v_equal = [];
end

grid_seed = [mu_grid, log(max(sig_grid, sig_min)), v_equal];
if length(grid_seed) == nParams
    seeds = [seeds; grid_seed];
end

% --- Strategy 3: Tight Cluster Split ---
sig_tight = repmat(std_global/2, 1, kComp);
mu_tight  = repmat(mu_global, 1, kComp) + randn(1,kComp)*std_global*0.1;
tight_seed = [sort(mu_tight), log(sig_tight), v_equal];

if length(tight_seed) == nParams
    seeds = [seeds; tight_seed];
end

% --- Safety Net ---
if isempty(seeds)
    mu_safe  = linspace(min(x), max(x), kComp);
    sig_safe = repmat(std_global, 1, kComp);
    if kComp > 1
        v_safe = zeros(1, kComp-1);
    else
        v_safe = [];
    end
    seeds = [mu_safe, log(sig_safe), v_safe];
end
end

%% -------------------- Parameter Reconstruction --------------------
function Params = ReconstructParams(p, kComp)
% Convert optimizer vector back to [mu, sigma, weight] matrix
mu = p(1:kComp);
sigma = exp(p(kComp+1:2*kComp));

if kComp > 1
    v = [p(2*kComp+1:end), 0];
    v = v - max(v); % Stable softmax
    w = exp(v) / sum(exp(v));
else
    w = 1;
end

% Sort by Mean for consistent output
Table = [mu(:), sigma(:), w(:)];
Params = sortrows(Table, 1);
end

%% -------------------- PDF Evaluation --------------------
function y = EvaluateGMM(x, Params)
y = zeros(size(x));
for k = 1:size(Params, 1)
    mu = Params(k,1);
    sig = Params(k,2);
    w = Params(k,3);
    y = y + w * normpdf(x, mu, sig);
end
end