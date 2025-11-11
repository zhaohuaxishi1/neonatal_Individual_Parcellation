%% analyze_correlation_variability_loading.m
% Correlate individual variability with network loading
% Author: jlzhao

clc; clear;

%% ------------------- Configuration ---------------------
NetNum = 11;

% Load group-level parcellation (initV)
load(fullfile('data', 'init.mat'), 'initV');

% Get hard network labels
[Network_Loading, network_labels] = max(initV, [], 2);  % Network_Loading: max loading values

% Load voxelwise individual variability
load(fullfile('output', 'individual_variability', 'voxelwise_variability.mat'), 'voxelwise_variability');
Individual_Variability = voxelwise_variability;

%% ------------------- Global Correlation ---------------------
[r_global, p_global] = corr(Network_Loading, Individual_Variability);

fprintf('Global correlation: r = %.4f, p = %.4f\n', r_global, p_global);


%% ------------------- Network-wise Correlation ---------------------
r_network = zeros(NetNum, 1);
p_network = zeros(NetNum, 1);

for i = 1:NetNum
    idx = find(network_labels == i);
    [r_network(i), p_network(i)] = corr(Network_Loading(idx), Individual_Variability(idx));
end

V_median = r_network;
V_median(:,2) = [11 5 10 3 4 9 8 7 6 2 1];

[V_median_sorted,V_median_sorted_index] = sort(V_median(:,1));
V_median_result = V_median(V_median_sorted_index,:);

%% ------------------- Save Results ---------------------
out_dir = fullfile('output', 'variability_loading_corr');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save(fullfile(out_dir, 'Individual_Variability.mat'), 'Individual_Variability');
save(fullfile(out_dir, 'Network_Loading.mat'), 'Network_Loading');
save(fullfile(out_dir, 'Correlation_Global.mat'), 'r_global', 'p_global');
save(fullfile(out_dir, 'Correlation_ByNetwork.mat'), 'r_network', 'p_network');


