function pi = stationary_dist(P)
%STATIONARY_DIST  Stationary distribution of a row-stochastic matrix P.
%   Solves the overdetermined system [P'-I; 1...1]*pi = [0...0; 1] in the
%   least-squares sense via QR (MATLAB \). More numerically stable than
%   eigendecomposition: no arbitrary fallback threshold, no iteration, and
%   well-conditioned even for near-reducible P (e.g. fly with no sleep
%   transitions). Negative numerical noise clamped; uniform fallback if
%   the solution collapses (degenerate input only).
%
%   Promoted to its own file on 2026-10-18. It had been duplicated as a local
%   function in four places. Three of those were mathematically identical for
%   any non-degenerate P and differed only in how they handled a collapsed
%   solution; this is the most defensive of them, from COMPUTE_HMM_ENERGETICS,
%   and it is the version the callers now share.
%
%   NOT adopted by IDENTIFYTRANSITIONDRIVERS, whose local copy divides by
%   max(sum(pi), 1e-12) instead of falling back to uniform. That is a real
%   behavioural difference in the degenerate case, so it was left alone rather
%   than changed silently.
    K  = size(P, 1);
    A  = [(P' - eye(K)); ones(1,K)];
    b  = [zeros(K,1);    1        ];
    pi = max(0, real(A \ b));
    s  = sum(pi);
    if s > realmin('double')
        pi = pi / s;
    else
        pi = ones(K,1) / K;   % uniform fallback for degenerate P only
    end
end
