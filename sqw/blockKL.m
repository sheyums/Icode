function D = blockKL(y, ph, nP, options)
%BLOCKKL  KL divergence between n-minute blocks and their time reverses, per phase (sqw_ref.block_kl; PRESPEC 20p).
%
%   D = blockKL(y, ph, nP)
%   D = blockKL(y, ph, nP, QW=false, NMax=4)
%
%   y    (T x 1) symbols 1 = Q, 2 = W, 3 = S, 0 = missing, 4 = unknown wake
%   ph   (T x 1) phase of each minute, 1..nP
%   QW   true: wake-only blocks (alphabet {Q, W}; a block containing S is not used)
%
%   D    (nP x NMax): D(p, n) = D_n = sum_b p(b) ln[p(b) / p(reverse b)], D_1 = 0.
%        p(b) = (c(b) + 0.5) / (N + 0.5 x alphabet^n): pseudocount 0.5 on every
%        block of the alphabet, so a block and its reverse that never occur get
%        equal mass and contribute 0. Blocks are n consecutive minutes, all in
%        phase p, none missing or unknown wake (overlapping blocks).
%
%   d_n = D_n - D_{n-1} (x 60 for nats/h); for a stationary first-order chain
%   d_3 = d_4 = d_2. Plug-in KL is biased upward with n: interpret only against
%   a matched null (BLOCKKLSTATS).
%
%   See also BLOCKKLSTATS.

arguments
    y (:,1) double
    ph (:,1) double {mustBePositive, mustBeInteger}
    nP (1,1) double {mustBePositive, mustBeInteger}
    options.QW (1,1) logical = false
    options.NMax (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.NMax, 1)} = 4
end
T = numel(y);
if numel(ph) ~= T, error('sqw:length', 'y and ph must have the same length.'); end
if options.QW, ok = y == 1 | y == 2; else, ok = y >= 1 & y <= 3; end
s = min(max(y - 1, 0), 2);                          % 0-based digit (only used where ok)
D = zeros(nP, options.NMax);
for n = 2:options.NMax
    m = T - n + 1;
    if m < 1, continue, end
    good = true(m, 1); code = zeros(m, 1);
    for k = 0:n - 1
        good = good & ok(1 + k:k + m) & ph(1 + k:k + m) == ph(1:m);
        code = code * 3 + s(1 + k:k + m);           % first minute most significant
    end
    allc = (0:3^n - 1)';
    digits = mod(floor(allc ./ 3.^(0:n - 1)), 3);   % digits(c+1, k+1) = coefficient of 3^k
    rev = digits(:, end:-1:1) * 3.^(0:n - 1)';      % code of the reversed block
    if options.QW, alph = all(digits < 2, 2); else, alph = true(3^n, 1); end
    for p = 1:nP
        sel = good & ph(1:m) == p;
        c = accumarray(code(sel) + 1, 1, [3^n 1]);
        pr = zeros(3^n, 1); pr(alph) = c(alph) + 0.5; pr = pr / sum(pr);
        D(p, n) = sum(pr(alph) .* log(pr(alph) ./ pr(rev(alph) + 1)));
    end
end
end
