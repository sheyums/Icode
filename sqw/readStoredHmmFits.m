function S = readStored20bFits(jsonlFile)
%READSTORED20BFITS  cClaude's stored 20b HMM fits (hmm_audit/sqw20b_fit2_theta0.5.jsonl) as MATLAB structs.
%
%   S = readStored20bFits(jsonlFile)
%
%   S  struct array, one per fly: g (genotype), fly0 (0-based index in sorted file
%      order, so the file is SQWLISTFILES(folder)(fly0 + 1)), T, deg, and the
%      LABELLED K = 3 fit as pi (3 x 1), P (3 x 3 x nP), B (3 x 3), state order
%      H1, H2, H3, symbols Q, W, S. Raw fits per K are in F (as jsondecode gives them).
%
%   See also CHECK20BEXACT, SQWMEMORYGROUP, SQWTURGROUP.

L = readlines(jsonlFile); L = L(strlength(L) > 0);
S = struct('g', {}, 'fly0', {}, 'T', {}, 'deg', {}, 'pi', {}, 'P', {}, 'B', {}, 'F', {});
for k = 1:numel(L)
    r = jsondecode(L(k));
    P = double(r.P3);
    if ismatrix(P) && isequal(size(P), [3 3]), nP = 1; else, nP = size(P, 1); P = permute(reshape(P, nP, 3, 3), [2 3 1]); end
    S(k) = struct('g', string(r.g), 'fly0', r.fly, 'T', r.T, 'deg', logical(r.deg), 'pi', double(r.pi3(:)), ...
        'P', P, 'B', double(reshape(r.B3, 3, 3)), 'F', r.F);
end
end
