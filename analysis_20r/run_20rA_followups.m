function Fu = run_20rA_followups(part)
% RUN_20RA_FOLLOWUPS  PRESPEC "2026-10-31 20r-A follow-ups". Run BY NAME from
% the Icode folder, with analysis_20r on the path.
%
%   part = 'halves' : (b) CS DD rise test within each 12-h half of the record clock
%   part = 'power'  : (a) 199 per0 DD subsamples at CS DD's n, rise test with
%                     NullReplicates 199; power per bin count

here = fileparts(mfilename('fullpath'));
switch part
    case 'halves'
        B = load(fullfile(here, 'CSDD_bouts_20r.mat'));
        on = B.metaW.OnsetBins;
        first = mod(on - 1, 86400) < 43200;
        W1 = B.Wake; W1(~first) = NaN;
        W2 = B.Wake; W2(first | isnan(on)) = NaN;
        Fu.first  = rise_test_20r(W1, B.xminWake, 'CS DD wake, record hours 0-12 (follow-up b)');
        Fu.second = rise_test_20r(W2, B.xminWake, 'CS DD wake, record hours 12-24 (follow-up b)');
        save(fullfile(here, 'followup_halves_CSDD.mat'), 'Fu', '-v7');

    case 'power'
        B = load(fullfile(here, 'per0DD_bouts_20r_replication.mat'));
        nCS = 3309; nRep = 199;
        keep = isfinite(B.Wake) & B.Wake >= B.xminWake;
        nFly = size(B.Wake, 1);
        bouts = arrayfun(@(i) B.Wake(i, keep(i, :)), (1:nFly)', 'UniformOutput', false);
        ni = cellfun(@numel, bouts);
        P = nan(nRep, 4); RR = nan(nRep, 4); verdict = strings(nRep, 1); nWhole = nan(nRep, 1);
        for b = 1:nRep
            rng(3000 + b, 'twister');
            ord = randperm(nFly);
            cs = cumsum(ni(ord));
            k = find(cs <= nCS, 1, 'last'); if isempty(k), k = 0; end
            x = [bouts{ord(1:k)}];
            rest = nCS - numel(x);
            if rest > 0
                nxt = bouts{ord(k + 1)};
                x = [x, nxt(randperm(numel(nxt), rest))]; %#ok<AGROW>
            end
            assert(numel(x) == nCS);
            nWhole(b) = k;
            R = rise_test_20r(x, B.xminWake, sprintf('per0 DD subsample %d at n = %d', b, nCS), 199);
            P(b, :) = R.RiseNullP; RR(b, :) = R.RiseRatio; verdict(b) = R.verdict;
        end
        power = mean(P < 0.05, 1);
        fprintf('\nFollow-up (a): power at n = %d, %d subsamples (whole flies per subsample: median %d)\n', ...
            nCS, nRep, median(nWhole));
        fprintf('  bins 12/16/20/24: power %.3f %.3f %.3f %.3f\n', power);
        fprintf('  median RiseRatio %.3f %.3f %.3f %.3f; median RiseNullP %.3f %.3f %.3f %.3f\n', ...
            median(RR), median(P));
        fprintf('  verdict "supported": %.3f; "supported" or "bin-sensitive": %.3f; NullDropped > 0: %d\n', ...
            mean(verdict == "supported"), ...
            mean(verdict == "supported" | verdict == "supported but bin-sensitive"), ...
            sum(startsWith(verdict, "not computable")));
        Fu = struct('nCS', nCS, 'P', P, 'RiseRatio', RR, 'verdict', verdict, 'nWhole', nWhole, 'power', power);
        save(fullfile(here, 'followup_power_per0DD.mat'), 'Fu', '-v7');
end
end
