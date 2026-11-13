% Export bout tables and MATLAB reference values for the Python reproduction kit.
% Every reference number is read from the saved results of the ledger entries
% (L, M, N, O, P, Q); nothing is recomputed here except group summaries of them.
sc  = 'C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad';
kit = fullfile(sc, 'bout_kit');
mk = @(p) mkdir(p);
for p = {fullfile(kit,'data','bouts'), fullfile(kit,'reference')}, if ~isfolder(p{1}), mk(p{1}); end, end

H  = load(fullfile(sc,'homeo_run','homeo_results.mat'), 'flies');
T0 = load(fullfile(sc,'N_run','tim0_flies.mat'), 'tim0');
D  = struct('per0DD', H.flies.per0DD, 'CSLD', H.flies.CSLD, 'CSDD', H.flies.CSDD);
D.tim0DD = T0.tim0;
names = fieldnames(D);

% ---------- 1. bout tables (whole records, edge bouts already removed) ----------
for k = 1:numel(names)
    F = D.(names{k}); rows = {};
    for f = 1:numel(F)
        n = numel(F(f).Dur);
        rows{end+1} = table(repmat(f, n, 1), repmat(string(F(f).Name), n, 1), ...
            double(F(f).Onset(:)), double(F(f).Dur(:)), double(F(f).IsSleep(:)), ...
            'VariableNames', {'fly_index','file','onset_s','dur_s','is_sleep'}); %#ok<SAGROW>
    end
    writetable(vertcat(rows{:}), fullfile(kit,'data','bouts',[names{k} '.csv']));
    fprintf('bouts %s: %d flies\n', names{k}, numel(F));
end

% ---------- 2. entry M / P / O(M): sleep-wake coupling ----------
Mp = load(fullfile(sc,'M_run','M_results.mat'), 'out');
Pp = load(fullfile(sc,'P_run','P_results.mat'), 'out');
Om = load(fullfile(sc,'O_run','OM_results.mat'), 'out');
MR = struct('per0DD', Mp.out.per0DD, 'CSLD', Mp.out.CSLD, 'CSDD', Pp.out.CSDD, 'tim0DD', Om.out.tim0DD);
cutH = struct('per0DD', 24, 'CSLD', 0, 'CSDD', 24, 'tim0DD', 24);
grp = struct();
for k = 1:numel(names)
    nm = names{k}; R = MR.(nm); F = D.(nm);
    t = array2table(R.obs, 'VariableNames', {'R','r_pre','r_post','dtail','reverse_tail','cpre','cpost'});
    t = addvars(t, string({F.Name}).', 'Before', 1, 'NewVariableNames', 'file');
    writetable(t, fullfile(kit,'reference',['coupling_perfly_' nm '.csv']));
    g = struct('dataset', nm, 'cut_hours', cutH.(nm), 'n_flies', numel(F), ...
        'nulls', {cellstr(R.nulls)}, 'permutations', size(R.Gn, 1), 'matlab_seed', 5, ...
        'R_observed', R.G(1), 'r_pre_observed', R.G(2), 'r_post_observed', R.G(3), ...
        'dtail_observed', R.G(4), 'reverse_tail_observed', R.G(5), ...
        'cpre_mean', mean(R.obs(:,6),'omitnan'), 'cpost_mean', mean(R.obs(:,7),'omitnan'));
    for j = 1:numel(R.nulls)
        g.(['R_null_mean_' char(R.nulls(j))]) = mean(R.Gn(:,1,j));
        g.(['R_null_sd_'   char(R.nulls(j))]) = std(R.Gn(:,1,j));
        g.(['R_p_'         char(R.nulls(j))]) = R.p(j);
        g.(['rpre_null_mean_' char(R.nulls(j))]) = mean(R.Gn(:,2,j));
    end
    grp.(nm) = g;
end
fid = fopen(fullfile(kit,'reference','coupling_group.json'), 'w'); fprintf(fid, '%s', jsonencode(grp, PrettyPrint=true)); fclose(fid);

% ---------- 3. entry L / O(L): long-bout clustering (per0 DD and tim0 DD, from 24 h) ----------
Lp = load(fullfile(sc,'L_run','L_results.mat'), 'R24');
Lt = load(fullfile(sc,'O_run','OL_results.mat'), 'R24');
LR = struct('per0DD', Lp.R24, 'tim0DD', Lt.R24); grpL = struct(); nul = ["P","S12","W2"];
for nm = ["per0DD","tim0DD"]
    R = LR.(nm); F = D.(nm);
    up = squeeze(R.obs(:,:,1)); lo = squeeze(R.obs(:,:,2));
    t = array2table([up lo], 'VariableNames', {'upper_lag1','upper_lag2','upper_lag3','upper_lag4', ...
        'lower_lag1','lower_lag2','lower_lag3','lower_lag4'});
    t = addvars(t, string({F.Name}).', 'Before', 1, 'NewVariableNames', 'file');
    writetable(t, fullfile(kit,'reference',['clustering_perfly_' char(nm) '.csv']));
    g = struct('dataset', nm, 'cut_hours', 24, 'upper_observed', R.G(1), 'lower_observed', R.G(2), ...
        'lag_profile_observed', R.gk, 'lag_profile_S12_null_mean', R.mk, 'Dlag', R.D, 'Dlag_p', R.pD);
    for j = 1:3
        g.(['upper_null_mean_' char(nul(j))]) = mean(R.Gn(:,1,j));
        g.(['upper_null_sd_' char(nul(j))])   = std(R.Gn(:,1,j));
        g.(['upper_p_' char(nul(j))])         = R.pt(1,j);
        g.(['lower_null_mean_' char(nul(j))]) = mean(R.Gn(:,2,j));
        g.(['lower_p_' char(nul(j))])         = R.pt(2,j);
    end
    grpL.(nm) = g;
end
fid = fopen(fullfile(kit,'reference','clustering_group.json'), 'w'); fprintf(fid, '%s', jsonencode(grpL, PrettyPrint=true)); fclose(fid);

% ---------- 4. entry Q: one-step link (pooled per0 DD + CS LD + CS DD) ----------
Q = load(fullfile(sc,'Q_run','Q_results.mat'));
ds = ["per0DD","CSLD","CSDD"]; files = [string({D.per0DD.Name}) string({D.CSLD.Name}) string({D.CSDD.Name})].';
t = array2table([Q.R.obs Q.R.nullmean Q.R.fast], 'VariableNames', ...
    {'r1','r2','r3','r_ss','r_ww','null_r1','null_r2','null_r3','null_ss','null_ww', ...
     'fast1','fast2','fast3','fast_ss','fast_ww'});
t = addvars(t, ds(Q.ds(:)).', files, 'Before', 1, 'NewVariableNames', {'dataset','file'});
writetable(t, fullfile(kit,'reference','onestep_perfly.csv'));
gq = struct('F', Q.R.F, 'kappa', Q.kap, 'kappa_ci', Q.ci, 'T1_p', Q.R.pT1, ...
    'mean_r2_observed', mean(Q.R.obs(isfinite(Q.R.obs(:,2)),2)), 'r2_null_group_mean', mean(Q.R.g2null), ...
    'r2_null_group_sd', std(Q.R.g2null), 'permutations', numel(Q.R.g2null), 'bootstraps', numel(Q.kb));
fid = fopen(fullfile(kit,'reference','onestep_group.json'), 'w'); fprintf(fid, '%s', jsonencode(gq, PrettyPrint=true)); fclose(fid);

% ---------- 5. entry N: first-day check (tim0 DD; per0 DD days 1-4) ----------
N = load(fullfile(sc,'N_run','N_results.mat'));
for k = 1:2
    o = N.out(k); nm = ["tim0DD","per0DD"]; F = D.(nm(k));
    L = reshape(permute(o.L, [1 3 2]), size(o.L,1), []); A = reshape(permute(o.A, [1 3 2]), size(o.A,1), []);
    vn = {}; for dd = 1:4, for hh = 1:2, vn{end+1} = sprintf('d%dh%d', dd, hh); end, end %#ok<SAGROW>
    t = [array2table(L, 'VariableNames', strcat('L_', vn)) array2table(A, 'VariableNames', strcat('A_', vn))];
    t = addvars(t, string({F.Name}).', 'Before', 1, 'NewVariableNames', 'file');
    writetable(t, fullfile(kit,'reference',['firstday_perfly_' char(nm(k)) '.csv']));
end
fprintf('reference written\n');
