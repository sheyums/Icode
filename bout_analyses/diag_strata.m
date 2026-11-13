cd('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts\Icode');
assert(isequal(fileparts(which('BoutIrreversibility')), pwd));
V = struct('name', {'A memory+profile, 1-h bins (as test 18)', 'B NO memory, profile, 1-h', ...
    'C memory, FLAT profile, 1-h', 'D memory+profile, 15-min bins', 'E memory+profile, 1-h, ZT positions FIXED (no timeline coupling)'}, ...
    'phi', {0.5, 0, 0.5, 0.5, 0.5}, 'flat', {false, false, true, false, false}, ...
    'binh', {1, 1, 1, 0.25, 1}, 'fixedzt', {false, false, false, false, true});
for v = V
  rej = [0 0]; nullMean = [];
  for i = 1:100
    xs = cell(1, 10); st = xs;
    for f = 1:10, [xs{f}, st{f}] = ldfly(v, 7000 + 100*i + f); end
    R = BoutIrreversibility(xs, Strata=st, Lags=1, NumSurrogates=199, MinPairs=10, Null="strata");
    rej = rej + ([R.Group.Bowley_P_Strata, R.Group.UpFraction_P_Strata] < 0.05);
    nullMean(end+1, :) = [R.Group.UpFraction, mean(R.GroupNull.strata.UpFraction)]; %#ok<AGROW>
  end
  fprintf('%-70s group rej Bowley %.2f Up %.2f | mean obs Up %.4f, mean null Up %.4f\n', v.name, rej/100, mean(nullMean));
end
function [sleep, strat] = ldfly(v, seed)
rs = RandStream('mt19937ar', 'Seed', seed); phi = v.phi;
t = 0; z = rs.randn; sleep = []; strat = []; k = 0;
while t < 5*86400
    k = k + 1;
    if v.fixedzt, tt = (k - 1) * 3200; else, tt = t; end   % position clock independent of durations
    zt = mod(tt/3600, 24);
    if v.flat, mw = 16; else, mw = 25*(zt < 12) + 8*(zt >= 12); end
    t = t + 2 + round(60*mw*-log(rs.rand));
    if ~v.fixedzt, tt = t; end
    zt = mod(tt/3600, 24);
    if v.flat, ms = 30; elseif zt < 12, ms = 15; else, ms = 20 + 60*(zt - 12)/12; end
    z = phi*z + sqrt(1 - phi^2)*rs.randn; e = -log(0.5*erfc(z/sqrt(2)));
    d = 300 + round(60*ms*e);
    sleep(end+1, 1) = d; strat(end+1, 1) = floor(zt / v.binh); %#ok<AGROW>
    t = t + d;
    if v.fixedzt && k*3200 > 5*86400, break, end
end
end
