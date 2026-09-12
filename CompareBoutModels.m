function R = CompareBoutModels(eventseries, xmin, options)
%COMPAREBOUTMODELS Fit a library of left-truncated models to bout durations
%and compare them by AICc and BIC, test the winner's absolute fit, and plot
%it against the data.
%
%   R = COMPAREBOUTMODELS(eventseries, xmin, SamplingInterval=dt, ...)
%
%   Pushes one dataset through every model in the library, ranks them,
%   checks whether the ranked winner actually FITS (ranking alone only says
%   which candidate is least bad), and walks down the ranking until one
%   passes. Returns the whole table, not just the winner.
%
%   THE LIBRARY
%     hyperexponential K=1..MaxComponents   FITHYPEREXPONENTIALMLE, one row
%                                           per K (see FLATTENING below)
%     exponentiated Weibull                 FITEXPONENTIATEDWEIBULLMLE
%     weibull                               FITWEIBULLMLE
%     gamma                                 FITGAMMAMLE
%     powerlaw_cutoff                       FITPOWERLAWCUTOFFMLE
%     pearson3                              FITPEARSON3MLE
%     powerlaw                              FITPOWERLAWMLE
%     erlang                                FITERLANGMLE
%     chisquared                            FITCHISQUAREDMLE
%     beta                                  FITBETAMLE, only if UpperBound
%                                           is supplied
%
%   FLATTENING. Every hyperexponential order K gets its own row and
%   competes directly against the other families, rather than the family
%   first picking a winner internally and only that winner being compared.
%   Two-stage selection hides the selection that already happened inside
%   the family; AICc and BIC compare MODELS, and a family boundary is not
%   special. Orders that fail the mixture fitter's identifiability gating
%   appear in the table with their reason and are barred from winning.
%
%   NESTING. Several library members are the same distribution, or special
%   cases of one another, and their rows are NOT independent evidence:
%     powerlaw_cutoff IS gamma, with alpha = 1-shape (identical logL, k,
%       AICc, BIC -- it is fitted once and relabelled);
%     hyperexponential K=1 IS the exponential, which is also gamma with
%       shape 1, weibull with shape 1, erlang with shape 1, and the
%       exponentiated Weibull with alpha=k=1;
%     weibull IS the exponentiated Weibull with alpha=1;
%     erlang and chisquared are both constrained gammas.
%   R.Nesting spells this out, and identical log-likelihoods among those
%   rows are a CHECK on the implementations rather than a coincidence: if
%   hyperexponential K=1 and a shape-1 gamma disagree by more than ~1e-9,
%   one of them is wrong.
%
%   COMPARABILITY. Every model is fitted with the same DistributionType,
%   the same SamplingInterval and the same xmin, on the same retained
%   observations, because information criteria may only be compared between
%   fits sharing a dominating measure. A discrete log-likelihood sits
%   roughly n*log(dt) above a continuous one, a difference of measure rather
%   than of evidence. This is enforced, not merely documented: the options
%   are passed to every fitter from one place.
%
%   AICc VERSUS BIC. Both are reported and they will often disagree, by
%   design: BIC's penalty is k*log(n), which at n=3400 is about 8.1 per
%   parameter against AIC's 2, so BIC systematically prefers fewer
%   components. Report both and say which you used. RankBy chooses which
%   orders the table and which winner is carried into the fit test; it does
%   not hide the other.
%
%   AKAIKE WEIGHTS. w_i = exp(-dAICc_i/2) / sum_j exp(-dAICc_j/2), the
%   relative support for each model given the candidate set. These matter
%   more than the winner's identity: if the top three split weight
%   0.4/0.35/0.25 there is no defensible single winner, and a bare
%   "model X won" would conceal that. Weights are computed over the
%   non-degenerate rows only.
%
%   GOODNESS OF FIT. Ranking is relative; a G test asks whether the model
%   is acceptable in absolute terms. Observations are binned (quantile
%   edges, then pooled until every expected count reaches
%   GoFMinExpected), and
%
%       G = 2 * sum_b O_b * log(O_b / E_b),   E_b = n*(Sh(lo_b) - Sh(hi_b))
%
%   using each fit's SurvivalHandle. With parameters estimated from
%   UNGROUPED data the null distribution of G is not exactly chi-square: it
%   lies between chi2(B-1-k) and chi2(B-1) (Chernoff & Lehmann 1954), so
%   both bounds are reported. GoFBootstrap="auto" additionally runs
%   a parametric bootstrap, which has the correct null distribution because
%   it refits the model to each replicate and so absorbs the estimation
%   effect that the df bounds merely bracket. It costs B refits, so "auto"
%   spends it only where it can change the answer: when the two bounds
%   STRADDLE GoFAlpha (pLower <= alpha < pUpper -- an outright undecided
%   test), or when pUpper lands in 0.01 to 0.2, near enough to a boundary
%   that the bracket's width matters. Elsewhere both bounds already agree.
%   Pass 0 to disable, or an integer B to force it; "auto" uses B=199.
%
%   The verdict rule follows from the same bracket. With no bootstrap,
%   Passed uses pLower, the conservative bound, so a model called acceptable
%   has passed the harder of the two readings. With a bootstrap, Passed uses
%   the bootstrap p. R.GoF.Ambiguous records whether the bounds disagreed.
%
%   Note that G depends on the binning. GoFBins and GoFMinExpected are
%   recorded in R.GoF so the number is reproducible rather than an artefact
%   of defaults, and the test has high power at large n -- on a few
%   thousand bouts a model can be rejected for a misfit too small to
%   matter. Read the magnitude, not only the p-value.
%
%   NAME-VALUE OPTIONS
%   SamplingInterval  REQUIRED, same units as xmin.
%   DistributionType  "discrete" (default) or "continuous", applied to all.
%   MaxComponents     highest hyperexponential order, default 4.
%   UpperBound        fixed upper support limit for beta, in DATA UNITS
%                     (86400 for a 24 h recording scored in seconds).
%                     Omit to leave beta out of the library.
%   Models            cellstr/string array to restrict the library.
%   RankBy            "AICc" (default) or "BIC".
%   GoFBins           number of quantile bins before pooling, default 40.
%   GoFMinExpected    pooling threshold on expected counts, default 5.
%   GoFAlpha          significance level for "passes", default 0.05.
%   GoFBootstrap      "auto" (default), 0, or a positive integer B.
%   Plot              logical, default true.
%   RandomSeed, Verbose
%
%   OUTPUT
%   R.Table        one row per candidate, ranked: Model, k, nReported, n,
%                  LogLik, AIC, AICc, BIC, dAICc, AkaikeWeight, dBIC,
%                  Degenerate, Reason, ParamText, ParamNames, Params,
%                  ParamSE. The last three vary in length between rows, so
%                  in MATLAB they are cell columns -- R.Table.Params{i} --
%                  while Octave, which has no table type, gets the
%                  underlying struct array and R.Table(i).Params.
%                  NOTE k and nReported differ wherever a
%                  constraint or an integer parameter is involved -- the
%                  hyperexponential at K=3 reports 6 numbers (3 tau, 3 q)
%                  but has k=5 free parameters because sum(q)=1 removes
%                  one. The information criteria use k; do not recompute
%                  them from numel(Params).
%   R.Fits         one record per row, parallel to R.Table: Model,
%                  SurvivalHandle, Refit (the reduced-start handle the
%                  bootstrap drives), and Full, which holds the originating
%                  fitter's own output untouched -- H.AllFits(K) for a
%                  hyperexponential row, the whole H for every other
%   R.Selected     the highest-ranked row that also passed the fit test,
%                  empty if none did
%   R.SelectedFit  that row's R.Fits record, or empty
%   R.SelectionPath what was skipped on the way down, and why
%   R.GoF          G, bins, df bounds, chi-square p at both, bootstrap p
%                  and B if run, Passed
%   R.Nesting      the equivalences listed above, as text
%   R.Figure       figure handle, or empty
%   R.n, R.xmin, R.SamplingInterval, R.DistributionType
%
%   See also FITHYPEREXPONENTIALMLE, FITEXPONENTIATEDWEIBULLMLE,
%   FITTRUNCATEDDISCRETEMLE, FITGAMMAMLE, FITPEARSON3MLE, FITBETAMLE.

arguments
    eventseries double {mustBeReal}
    xmin (1,1) double {mustBePositive}
    options.SamplingInterval (1,1) double = NaN
    options.DistributionType (1,1) string {mustBeMember(options.DistributionType,["continuous","discrete"])} = "discrete"
    options.MaxComponents (1,1) double {mustBeInteger,mustBePositive} = 4
    options.UpperBound (1,1) double = NaN
    options.Models = []
    options.RankBy (1,1) string {mustBeMember(options.RankBy,["AICc","BIC"])} = "AICc"
    options.GoFBins (1,1) double {mustBeInteger,mustBePositive} = 40
    options.GoFMinExpected (1,1) double {mustBePositive} = 5
    options.GoFAlpha (1,1) double {mustBePositive} = 0.05
    options.GoFBootstrap = "auto"
    options.Plot (1,1) logical = true
    options.RandomSeed = []
    options.Verbose (1,1) logical = true
end

if isnan(options.SamplingInterval)
    error('CompareBoutModels:SamplingIntervalRequired', ...
        ['SamplingInterval is required, in the SAME UNITS as xmin. Every ' ...
         'model is fitted with it so the information criteria stay ' ...
         'comparable.']);
end
if ~isempty(options.RandomSeed)
    rng(options.RandomSeed);
end

dt = options.SamplingInterval;
mode = options.DistributionType;
common = {'SamplingInterval', dt, 'DistributionType', mode, 'Verbose', false};

% ------------------------------------------------------------- the library
cands = buildLibrary(options, common);
if ~isempty(options.Models)
    want = asCellstr(options.Models);
    cands = cands(ismember({cands.Name}, want));
    if isempty(cands)
        error('CompareBoutModels:NoModelsSelected', ...
            'None of the requested models are in the library.');
    end
end

if options.Verbose
    fprintf('CompareBoutModels: %d candidate(s), %s mode, xmin=%.6g, dt=%.6g.\n', ...
        numel(cands), mode, xmin, dt);
end

% ------------------------------------------------------------------ fitting
rows = struct([]); fits = struct([]); nr = 0;
for ci = 1:numel(cands)
    c = cands(ci);
    try
        out = c.Fit(eventseries, xmin);
    catch err
        if options.Verbose
            fprintf('  %-22s skipped: %s\n', c.Name, err.message);
        end
        continue
    end
    for oi = 1:numel(out)
        nr = nr + 1;
        rows(nr) = out(oi).Row;
        fits(nr) = out(oi).Fit;
    end
end
if nr == 0
    error('CompareBoutModels:NoFits', 'No model in the library could be fitted.');
end
R.n = rows(1).n;

% ------------------------------------------------- ranking and the weights
ok = ~[rows.Degenerate] & isfinite([rows.AICc]) & isfinite([rows.BIC]);
if ~any(ok)
    error('CompareBoutModels:AllDegenerate', ...
        ['Every candidate was flagged degenerate or produced a ' ...
         'non-finite criterion. Inspect R.Table.Reason.']);
end
aicc = [rows.AICc]; bicv = [rows.BIC];
bestA = min(aicc(ok)); bestB = min(bicv(ok));
wraw = zeros(1, nr);
wraw(ok) = exp(-(aicc(ok) - bestA)/2);
for i = 1:nr
    rows(i).dAICc = aicc(i) - bestA;
    rows(i).dBIC  = bicv(i) - bestB;
    rows(i).AkaikeWeight = wraw(i) / sum(wraw);
end

key = [rows.AICc];
if strcmp(options.RankBy, "BIC"), key = [rows.BIC]; end
key(~ok) = Inf;                       % degenerate rows rank last
[~, order] = sort(key);
rows = rows(order); fits = fits(order); ok = ok(order);

R.Table = rowsToTable(rows);
R.Fits = fits;
R.Nesting = nestingNotes();
R.xmin = xmin; R.SamplingInterval = dt; R.DistributionType = mode;
R.RankBy = options.RankBy;

% -------------------------------------------- goodness of fit, walking down
data = retainedData(eventseries, xmin, dt, mode);
path = {}; sel = []; gof = [];
for i = 1:nr
    if ~ok(i)
        path{end+1} = sprintf('%s: skipped, %s', rows(i).Model, rows(i).Reason); %#ok<AGROW>
        continue
    end
    g = gTest(data, xmin, dt, fits(i), rows(i).k, options);
    if g.Passed
        sel = i; gof = g;
        path{end+1} = sprintf('%s: ranked %d, G=%.2f p=%.4g -- accepted', ...
            rows(i).Model, i, g.G, g.pUpper); %#ok<AGROW>
        break
    end
    if isnan(g.G)
        path{end+1} = sprintf('%s: ranked %d, %s', ...
            rows(i).Model, i, g.Basis); %#ok<AGROW>
    else
        path{end+1} = sprintf('%s: ranked %d, G=%.2f p=%.4g -- rejected at alpha=%g', ...
            rows(i).Model, i, g.G, g.pUpper, options.GoFAlpha); %#ok<AGROW>
    end
    if isempty(gof), gof = g; end
end
R.SelectionPath = path(:);
R.GoF = gof;
if isempty(sel)
    R.Selected = [];
    R.SelectedFit = [];
    if options.Verbose
        fprintf(['CompareBoutModels: no candidate passed the fit test at ' ...
            'alpha=%g. R.Selected is empty; see R.SelectionPath.\n'], options.GoFAlpha);
    end
else
    R.Selected = rows(sel);
    R.SelectedFit = fits(sel);
end

% ----------------------------------------------------------------- reporting
if options.Verbose
    fprintf('\n%-22s %4s %12s %11s %11s %8s %7s\n', ...
        'model','k','logL','AICc','BIC','dAICc','w');
    for i = 1:nr
        flag = ''; if ~ok(i), flag = '  [excluded]'; end
        fprintf('%-22s %4d %12.2f %11.2f %11.2f %8.2f %7.3f%s\n', ...
            rows(i).Model, rows(i).k, rows(i).LogLik, rows(i).AICc, ...
            rows(i).BIC, rows(i).dAICc, rows(i).AkaikeWeight, flag);
    end
    if ~isempty(sel)
        fprintf('\nselected: %s   %s\n', rows(sel).Model, rows(sel).ParamText);
        fprintf('fit test: G=%.3f, %d bins, chi2 p in [%.4g, %.4g]', ...
            gof.G, gof.nBins, gof.pLower, gof.pUpper);
        if ~isnan(gof.pBootstrap)
            fprintf(', bootstrap p=%.4g (B=%d)', gof.pBootstrap, gof.B);
        end
        fprintf('\n');
    end
end

% -------------------------------------------------------------------- plot
R.Figure = [];
if options.Plot && ~isempty(sel)
    try
        R.Figure = plotFit(data, xmin, dt, fits(sel), rows(sel), gof, options);
    catch err
        if options.Verbose
            fprintf('CompareBoutModels: plotting skipped (%s).\n', err.message);
        end
    end
end
end

% ========================================================================
%                             THE LIBRARY
% ========================================================================
function cands = buildLibrary(options, common)
% One entry per FAMILY. Each .Fit call returns an array of {Row, Fit}
% records, so a family may contribute several rows -- the hyperexponential
% contributes one per order K. Models= filters at this family level, so
% "hyperexponential" selects every K rather than one of them.

names = {}; fitters = {};

mc = options.MaxComponents;
names{end+1}   = 'hyperexponential';
fitters{end+1} = @(d,x) fitHyperCand(d, x, mc, common);

names{end+1}   = 'exp_weibull';
fitters{end+1} = @(d,x) fitEWCand(d, x, common);

% Engine-backed families. Column 3 holds extra options for that family.
simple = { ...
    'weibull',         @FitWeibullMLE,        {}; ...
    'gamma',           @FitGammaMLE,          {}; ...
    'powerlaw_cutoff', @FitPowerLawCutoffMLE, {}; ...
    'pearson3',        @FitPearson3MLE,       {}; ...
    'powerlaw',        @FitPowerLawMLE,       {}; ...
    'erlang',          @FitErlangMLE,         {}; ...
    'chisquared',      @FitChiSquaredMLE,     {}};
if ~isnan(options.UpperBound)
    % Beta needs a finite upper support limit and has no sensible default:
    % the recording length is a choice about the experiment, not about the
    % data, so it is supplied or beta stays out of the library entirely.
    simple(end+1,:) = {'beta', @FitBetaMLE, {'UpperBound', options.UpperBound}};
end
for i = 1:size(simple,1)
    nm = simple{i,1}; fh = simple{i,2}; ex = [common, simple{i,3}];
    names{end+1}   = nm;                                        %#ok<AGROW>
    fitters{end+1} = @(d,x) fitSimpleCand(d, x, nm, fh, ex);     %#ok<AGROW>
end

cands = struct('Name', names, 'Fit', fitters);
end

% ------------------------------------------------------------------------
function out = fitSimpleCand(d, x, name, fitter, extra)
% Engine-backed families all report ParamNames/Params/ParamSE already, and
% for them every parameter is free, so nReported == k.
H = fitter(d, x, extra{:}, 'ErrorOnNoValidFit', false);
pn = asCellstr(H.ParamNames);
pv = H.Params(:).';
pse = H.ParamSE(:).';
row = makeRow(name, H.k, numel(pv), H.n, H.LogLik, H.AIC, H.AICc, H.BIC, ...
    H.Failed, H.FailureReason, pn, pv, pse);
rec = makeRec(name, H.SurvivalHandle, H, @(nd) refitSimple(nd, x, fitter, extra));
out = struct('Row', row, 'Fit', rec);
end

function r = refitSimple(nd, x, fitter, extra)
H = fitter(nd, x, extra{:}, 'ErrorOnNoValidFit', false, 'nStarts', 4);
r = struct('LogLik', H.LogLik, 'SurvivalHandle', H.SurvivalHandle, ...
    'Failed', H.Failed);
end

% ------------------------------------------------------------------------
function out = fitEWCand(d, x, common)
H = FitExponentiatedWeibullMLE(d, x, common{:}, 'ErrorOnNoValidFit', false);
pn  = {'lambda','k','alpha'};
pv  = [H.Lambda, H.K, H.Alpha];
pse = [H.LambdaSE, H.KSE, H.AlphaSE];
% FixAlpha would drop k to 2 while still reporting three numbers (alpha
% being pinned at 1, not estimated); nReported stays 3 either way.
row = makeRow('exp_weibull', H.k, 3, H.n, H.LogLik, H.AIC, H.AICc, ...
    H.k*log(H.n) - 2*H.LogLik, H.Failed, H.FailureReason, pn, pv, pse);
rec = makeRec('exp_weibull', H.SurvivalHandle, H, @(nd) refitEW(nd, x, common));
out = struct('Row', row, 'Fit', rec);
end

function r = refitEW(nd, x, common)
H = FitExponentiatedWeibullMLE(nd, x, common{:}, 'ErrorOnNoValidFit', false, ...
    'nStartsBase', 4, 'nStartsPerParameter', 2, 'maxStarts', 10);
r = struct('LogLik', H.LogLik, 'SurvivalHandle', H.SurvivalHandle, ...
    'Failed', H.Failed);
end

% ------------------------------------------------------------------------
function out = fitHyperCand(d, x, mc, common)
% One row per order. ErrorOnNoValidFit=false so that a family in which no
% order is identifiable still returns its rows, with reasons, instead of
% aborting the whole comparison.
H = FitHyperexponentialMLE(d, x, common{:}, 'MaxComponents', mc, ...
    'ErrorOnNoValidFit', false);
out = struct('Row', {}, 'Fit', {});
for K = 1:numel(H.AllFits)
    f = H.AllFits(K);
    name = sprintf('hyperexp K=%d', K);
    kfree = 2*K - 1;        % K rates + K weights, less sum(q)=1
    nrep  = 2*K;            % ...but 2K numbers are worth reporting
    pn = hyperNames(K);
    if ~f.Success
        why = f.DegenerateReason;
        if isempty(why), why = 'not fitted'; end
        row = makeRow(name, kfree, nrep, H.n, NaN, NaN, NaN, NaN, ...
            true, why, pn, nan(1,nrep), nan(1,nrep));
        rec = makeRec(name, [], f, []);
    else
        pv  = [f.Tau(:).',   f.WeightsObserved(:).'];
        pse = [f.TauSE(:).', f.WeightsObservedSE(:).'];
        row = makeRow(name, kfree, nrep, f.n, f.LogLik, f.AIC, f.AICc, ...
            kfree*log(f.n) - 2*f.LogLik, f.Degenerate, f.DegenerateReason, ...
            pn, pv, pse);
        rec = makeRec(name, f.SurvivalHandle, f, @(nd) refitHyper(nd, x, K, common));
    end
    out(end+1,1) = struct('Row', row, 'Fit', rec);                 %#ok<AGROW>
end
end

function r = refitHyper(nd, x, K, common)
% MaxComponents=K fits 1..K and selects; for the bootstrap we want order K
% specifically, so read AllFits(K) rather than H.Selected.
H = FitHyperexponentialMLE(nd, x, common{:}, 'MaxComponents', K, ...
    'ErrorOnNoValidFit', false, 'nStartsBase', 3, ...
    'nStartsPerComponent', 3, 'maxStarts', 12);
f = H.AllFits(K);
r = struct('LogLik', f.LogLik, 'SurvivalHandle', f.SurvivalHandle, ...
    'Failed', ~f.Success);
end

function pn = hyperNames(K)
pn = cell(1, 2*K);
for j = 1:K, pn{j}     = sprintf('tau%d', j); end
for j = 1:K, pn{K+j}   = sprintf('q%d',   j); end
end

% ========================================================================
%                          ROWS AND THE TABLE
% ========================================================================
function row = makeRow(model, k, nrep, n, logL, aic, aicc, bic, degen, reason, pn, pv, pse)
% Fixed field set and order, so rows from different families concatenate.
% Params carries ALL parameters, including ones a constraint makes
% redundant: at K=3 the hyperexponential reports three taus and all three
% q, even though sum(q)=1 means only two q are free. k stays 5 and the
% information criteria use k, never numel(Params) -- see the header.
if isempty(reason), reason = ''; end
if isempty(pn), pn = {}; end
row = struct( ...
    'Model',        model, ...
    'k',            k, ...
    'nReported',    nrep, ...
    'n',            n, ...
    'LogLik',       logL, ...
    'AIC',          aic, ...
    'AICc',         aicc, ...
    'BIC',          bic, ...
    'dAICc',        NaN, ...
    'AkaikeWeight', NaN, ...
    'dBIC',         NaN, ...
    'Degenerate',   logical(degen), ...
    'Reason',       reason, ...
    'ParamText',    paramText(model, pn, pv));
% Assigned rather than passed to struct(), which would read a cell value as
% a request for a struct ARRAY. Stored unwrapped: a struct array's fields
% may differ in length between elements, and only the table conversion
% needs them boxed.
row.ParamNames = pn(:).';
row.Params     = pv(:).';
row.ParamSE    = pse(:).';
end

function rec = makeRec(name, Sh, full, refit)
rec = struct('Model', name, 'SurvivalHandle', {Sh}, 'Refit', {refit}, ...
    'Full', {full});
end

function s = paramText(model, pn, pv)
if isempty(pv) || ~any(isfinite(pv))
    s = '(no fit)'; return
end
if strncmp(model, 'hyperexp', 8) && mod(numel(pv), 2) == 0
    K = numel(pv)/2;
    s = sprintf('tau=[%s], q=[%s]', ...
        strtrim(sprintf('%.4g ', pv(1:K))), ...
        strtrim(sprintf('%.4g ', pv(K+1:end))));
    return
end
parts = cell(1, numel(pv));
for i = 1:numel(pv)
    if i <= numel(pn), nm = pn{i}; else, nm = sprintf('p%d', i); end
    parts{i} = sprintf('%s=%.4g', nm, pv(i));
end
s = strjoin(parts, ', ');
end

function c = asCellstr(x)
% Accepts a cellstr, a char row, or a string array, without depending on
% the string class -- Octave has none, and the test harness runs there.
if iscell(x)
    c = cell(1, numel(x));
    for i = 1:numel(x), c{i} = char(x{i}); end
elseif ischar(x)
    c = {x};
else
    c = cellstr(x);
end
c = c(:).';
end

function tf = isAuto(v)
if isnumeric(v) || islogical(v), tf = false; return; end
tf = strcmpi(char(v), 'auto');
end

function T = rowsToTable(rows)
% struct2table needs same-height fields, so the variable-length parameter
% vectors travel as 1x1 cells (makeRow already wrapped them). Without the
% toolbox-free table type -- Octave -- hand back the struct array, which
% every field access in the header still works on.
if ~(exist('struct2table', 'file') == 2 || exist('struct2table', 'builtin') == 5)
    T = rows; return
end
% struct2table needs every field the same height, so the three
% variable-length parameter fields are boxed into cell columns here and
% only here. table2struct unboxes them again, which is what makes
% R.Table(i).Params and the Octave struct array's field agree in shape.
plain = rmfield(rows, {'ParamNames', 'Params', 'ParamSE'});
T = struct2table(plain);
T.ParamNames = {rows.ParamNames}.';
T.Params     = {rows.Params}.';
T.ParamSE    = {rows.ParamSE}.';
end

function notes = nestingNotes()
notes = { ...
 'powerlaw_cutoff IS gamma with alpha = 1 - shape: same fit, logL, k, AICc, BIC.'; ...
 'hyperexp K=1 IS the exponential = gamma(shape 1) = weibull(shape 1) = erlang(shape 1) = exp_weibull(alpha=k=1).'; ...
 'weibull IS exp_weibull with alpha = 1.'; ...
 'erlang and chisquared are both gammas under a constraint (integer shape; scale 2, shape nu/2).'; ...
 'Rows above are NOT independent evidence. Equal log-likelihoods among them are a check on the implementations, not a coincidence: a disagreement beyond ~1e-9 means one is wrong.'};
end

% ========================================================================
%                        DATA AND GOODNESS OF FIT
% ========================================================================
function d = retainedData(raw, xmin, dt, mode)
% Must reproduce the fitters' own retention rule exactly, or the expected
% counts are computed against a different n than the likelihood used. In
% discrete mode the support is round(t/dt) >= n_min, so filtering on raw
% t >= xmin would wrongly drop observations in [(n_min-1/2)*dt, xmin).
raw = raw(:);
raw = raw(isfinite(raw) & raw > 0);
if strcmp(mode, "discrete")
    n_min = max(1, round(xmin / dt));
    d = raw(round(raw / dt) >= n_min);
else
    d = raw(raw >= xmin);
end
end

% ------------------------------------------------------------------------
function g = gTest(data, xmin, dt, rec, k, options)
isD = strcmp(options.DistributionType, "discrete");
Sh  = rec.SurvivalHandle;

[O, E, edges, n, sumP] = binAndExpect(data, xmin, dt, Sh, isD, options);
B = numel(O);

% A fit whose survival handle puts no mass on the observed range is not
% something to compute a G statistic for. Pooling would collapse every bin
% into one, B-1 would be 0, and G would come out as a tidy 0.00 with a NaN
% p-value -- a broken fit looking like a perfect one. Say so instead.
if ~(sumP > 1 - 1e-6) || B < 2
    g = struct('G', NaN, 'nBins', B, 'dfLower', NaN, 'dfUpper', NaN, ...
        'pLower', NaN, 'pUpper', NaN, 'Ambiguous', false, ...
        'pBootstrap', NaN, 'B', 0, 'Passed', false, ...
        'Basis', sprintf(['not testable: the fitted survival puts %.4g of ' ...
            'its probability on the observed range, over %d bin(s)'], sumP, B), ...
        'Alpha', options.GoFAlpha, 'Observed', O, 'Expected', E, ...
        'Edges', edges(:).', 'ProbabilityMass', sumP, ...
        'GoFBins', options.GoFBins, 'GoFMinExpected', options.GoFMinExpected);
    return
end

G = gstat(O, E);

dfUpper = B - 1;                       % Chernoff & Lehmann (1954) bounds
dfLower = B - 1 - k;
pUpper = chi2Upper(G, dfUpper);        % the larger, laxer p
if dfLower > 0
    pLower = chi2Upper(G, dfLower);    % the smaller, conservative p
else
    pLower = NaN;
end

alpha = options.GoFAlpha;
ambiguous = isfinite(pLower) && pLower <= alpha && pUpper > alpha;

% --- does the bootstrap earn its keep here?
Bboot = 0;
if isnumeric(options.GoFBootstrap)
    Bboot = round(options.GoFBootstrap);
elseif isAuto(options.GoFBootstrap)
    if ambiguous || (pUpper >= 0.01 && pUpper <= 0.2)
        Bboot = 199;
    end
end
if isempty(rec.Refit), Bboot = 0; end   % nothing to refit with

pBoot = NaN; nValid = 0;
if Bboot > 0
    if options.Verbose
        fprintf('  bootstrap: %s, B=%d ...', rec.Model, Bboot);
    end
    Gb = nan(1, Bboot);
    for b = 1:Bboot
        try
            sim = sampleFromSurvival(Sh, n, xmin, dt, isD);
            rb  = rec.Refit(sim);
            if rb.Failed || isempty(rb.SurvivalHandle), continue; end
            [Ob, Eb] = binAndExpect(sim, xmin, dt, rb.SurvivalHandle, isD, options);
            Gb(b) = gstat(Ob, Eb);
        catch
            % a replicate that will not refit is dropped, not counted
        end
    end
    valid = Gb(isfinite(Gb));
    nValid = numel(valid);
    if nValid >= 20
        % +1 in both places: the observed statistic is itself one draw
        % from the null, which keeps the p-value from ever being 0.
        pBoot = (1 + nnz(valid >= G)) / (nValid + 1);
    end
    if options.Verbose, fprintf(' p=%.4g (%d valid)\n', pBoot, nValid); end
end

if isfinite(pBoot)
    passed = pBoot > alpha;
    basis  = 'bootstrap';
elseif isfinite(pLower)
    passed = pLower > alpha;
    basis  = 'chi2 at df=B-1-k (conservative bound)';
else
    passed = pUpper > alpha;
    basis  = 'chi2 at df=B-1 (too few bins for the lower bound)';
end

g = struct('G', G, 'nBins', B, 'dfLower', dfLower, 'dfUpper', dfUpper, ...
    'pLower', pLower, 'pUpper', pUpper, 'Ambiguous', ambiguous, ...
    'pBootstrap', pBoot, 'B', nValid, 'Passed', passed, 'Basis', basis, ...
    'Alpha', alpha, 'Observed', O, 'Expected', E, 'Edges', edges(:).', ...
    'ProbabilityMass', sumP, 'GoFBins', options.GoFBins, ...
    'GoFMinExpected', options.GoFMinExpected);
end

function G = gstat(O, E)
pos = O > 0;
G = 2 * sum(O(pos) .* log(O(pos) ./ E(pos)));
end

function p = chi2Upper(x, d)
% P(chi2_d > x) without the Statistics toolbox: the regularized upper
% incomplete gamma, which both MATLAB and Octave provide.
if ~(d > 0) || ~isfinite(x) || x < 0
    p = NaN; return
end
p = gammainc(x/2, d/2, 'upper');
end

% ------------------------------------------------------------------------
function [O, E, edges, n, sumP] = binAndExpect(data, xmin, dt, Sh, isD, options)
% Quantile bins, then pool adjacent bins until every EXPECTED count clears
% GoFMinExpected. Pooling on expected (not observed) counts is what the
% chi-square approximation actually requires.
n = numel(data);
if isD
    xs = round(data(:) / dt);           % work on the integer grid
    lo0 = max(1, round(xmin / dt));
else
    xs = data(:);
    lo0 = xmin;
end
xsort = sort(xs);

nb = min(options.GoFBins, max(2, floor(n / max(1, options.GoFMinExpected))));
qs = (1:nb-1) / nb;
eint = xsort(max(1, min(n, round(qs * n))));
edges = unique([lo0; eint(:); Inf]);
if numel(edges) < 3
    edges = [lo0; Inf];                 % degenerate data: one bin
end

% Observed. Bin b is [edges(b), edges(b+1)) -- on the integer grid that is
% the integers edges(b) .. edges(b+1)-1.
B = numel(edges) - 1;
O = zeros(1, B);
for b = 1:B
    O(b) = nnz(xs >= edges(b) & xs < edges(b+1));
end

% Expected. Both survival handles are P(T > t | T >= xmin) evaluated on the
% grid, so the cut for integer edge e is (e-1)*dt: S there is P(N >= e).
E = n * binProb(edges, Sh, isD, dt);

% --- pool
b = 1;
while numel(E) > 1
    if E(b) >= options.GoFMinExpected
        b = b + 1;
        if b > numel(E), break; end
        continue
    end
    if b == numel(E), lo = b - 1; else, lo = b; end
    E(lo) = E(lo) + E(lo+1);
    O(lo) = O(lo) + O(lo+1);
    E(lo+1) = []; O(lo+1) = []; edges(lo+1) = [];
    b = lo;
end
sumP = sum(E) / n;      % should be 1; a shortfall means the handle is off
end

function p = binProb(edges, Sh, isD, dt)
if isD
    cuts = (edges - 1) * dt;
else
    cuts = edges;
end
fin = isfinite(cuts);
Sv = zeros(numel(cuts), 1);
Sv(fin) = Sh(cuts(fin));
Sv(~fin) = 0;                  % S(Inf) = 0, and Sh may not accept Inf
p = max(Sv(1:end-1) - Sv(2:end), 0).';
end

% ------------------------------------------------------------------------
function s = sampleFromSurvival(Sh, n, xmin, dt, isD)
% Inverse-CDF draw from the FITTED truncated model using only its survival
% handle, so one sampler serves all three fitters. Both branches bisect
% rather than tabulate: enumerating the grid until the survival is
% negligible is unbounded work on a heavy tail (a fitted power law with
% alpha 0.8 needs ~1e15 grid steps to reach 1e-12), while bisection costs
% about 52 vectorized evaluations whatever the tail does.
u = rand(n, 1);
u(u <= 0) = eps;                        % tgt < 1 keeps the bracket valid
tgt = 1 - u;                            % want the quantile where Sh == tgt

if isD
    % Smallest integer j with Sh(j*dt) <= tgt, i.e. F(j) >= u. Sh is a
    % right-continuous step function on the grid, so the answer is exact.
    nmin = max(1, round(xmin / dt));
    hiS = nmin + 1;
    while Sh(hiS*dt) > min(tgt) && hiS < 2^52
        hiS = 2 * hiS;
    end
    lo = repmat(nmin - 1, n, 1);        % Sh(lo*dt) == 1 > tgt
    hg = repmat(hiS, n, 1);             % Sh(hg*dt) <= tgt
    while any(hg - lo > 1)
        mid = floor((lo + hg) / 2);
        above = Sh(mid * dt) > tgt;
        lo(above)  = mid(above);
        hg(~above) = mid(~above);
    end
    s = hg * dt;
else
    hi = 2 * xmin;
    while Sh(hi) > min(tgt) && hi < xmin * 1e12
        hi = 2 * hi;
    end
    lo = repmat(xmin, n, 1); hg = repmat(hi, n, 1);
    for it = 1:80
        mid = 0.5 * (lo + hg);
        above = Sh(mid) > tgt;          % Sh decreasing: t still too small
        lo(above)  = mid(above);
        hg(~above) = mid(~above);
    end
    s = 0.5 * (lo + hg);
end
end

% ========================================================================
%                                 PLOT
% ========================================================================
function fh = plotFit(data, xmin, dt, rec, row, gof, options)
isD = strcmp(options.DistributionType, "discrete");
n = numel(data);
xs = sort(data(:));

% Empirical survival P(T > t). No censoring here, so this is just 1-ECDF;
% a Kaplan-Meier estimator would reduce to the same thing.
ux = unique(xs);
cnt = zeros(numel(ux), 1);
for i = 1:numel(ux), cnt(i) = nnz(xs == ux(i)); end
Semp = 1 - cumsum(cnt) / n;

fh = figure('Name', sprintf('CompareBoutModels: %s', row.Model));

% --- survival, log-log: where a heavy tail either is or is not straight
subplot(2, 1, 1);
tg = logspace(log10(xmin), log10(max(xs) * 1.5), 400).';
if isD
    tg = unique(round(tg / dt) * dt);
    tg = tg(tg >= (max(1, round(xmin/dt)) - 1) * dt);
end
Sfit = rec.SurvivalHandle(tg);
kp = Semp > 0;
stairs([xmin; ux(kp)], [1; Semp(kp)], 'k-', 'LineWidth', 1.0); hold on
plot(tg, Sfit, 'r-', 'LineWidth', 1.6);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('bout duration'); ylabel('P(T > t | T \geq xmin)');
title(sprintf('%s   (n=%d, %s, dt=%g)   %s', row.Model, n, ...
    options.DistributionType, dt, row.ParamText), 'Interpreter', 'none');
legend({'data', 'fit'}, 'Location', 'southwest'); legend boxoff
grid on

% --- Pearson residuals from the very bins the G test used, so the plot and
% the p-value cannot disagree about where the misfit is
subplot(2, 1, 2);
r = (gof.Observed - gof.Expected) ./ sqrt(gof.Expected);
bar(1:numel(r), r, 'FaceColor', [0.35 0.35 0.7], 'EdgeColor', 'none');
hold on
plot(xlim, [ 2  2], 'k:'); plot(xlim, [-2 -2], 'k:');
xlabel(sprintf('G-test bin (%d bins, edges in R.GoF.Edges)', gof.nBins));
ylabel('(O-E)/sqrt(E)');
if isfinite(gof.pBootstrap)
    title(sprintf('G=%.2f, bootstrap p=%.4g (B=%d)', gof.G, gof.pBootstrap, gof.B));
else
    title(sprintf('G=%.2f, chi2 p in [%.4g, %.4g]', gof.G, gof.pLower, gof.pUpper));
end
grid on
end
