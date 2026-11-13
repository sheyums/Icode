L = load('..\homeo_run\homeo_results.mat', 'flies'); F = L.flies.CSLD;
for bw = [1 2 3]
    nb = []; for f = 1:numel(F)
        on = F(f).Onset(F(f).IsSleep); zt = mod(on,86400)/3600; k = zt < 12; on = on(k); zt = zt(k);
        g = floor(on/86400)*100 + floor(zt/bw); [~,~,gi] = unique(g); c = accumarray(gi,1); nb = [nb; c(gi)]; %#ok<AGROW>
    end
    fprintf('CS LD day, day x %d-h ZT bins: bouts in bins of size 1: %.2f, >=3: %.2f, median bin size per bout %g\n', bw, mean(nb==1), mean(nb>=3), median(nb));
end
