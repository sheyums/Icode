cd('C:\Users\sheyu\Box\MATLAB scripts\Rockefeller Matlab scripts');
fprintf('EstimatePhaseStamp: %s\n', which('EstimatePhaseStamp'));
D = 'C:\Users\sheyu\Box\HMM paper_Activity based\different genotypes activity data\per0 ld';
[T, figs] = EstimatePhaseStamp(D, ClockState="arrhythmic");
disp(T(:, setdiff(T.Properties.VariableNames, {'Files','Folder'}, 'stable')))
for k = 1:numel(figs), try exportgraphics(figs(k), fullfile('C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad', sprintf('eps_per0ld_%d.png', k)), Resolution=110); catch, end, end
save('C:\Users\sheyu\AppData\Local\Temp\claude\6d64e3807eb78b57\c516ceb9-b449-47ae-9473-3050e269184c\scratchpad\eps_per0ld.mat', 'T');
