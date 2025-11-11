%% plot_individual_variability_by_network.m
% Visualize network-level individual variability using violin plots
% Author: jlzhao

clc; clear;

%% ------------------- Setup ---------------------
addpath(genpath('utils/violinplot'));

% Load voxelwise variability scores and network assignment
load(fullfile('output', 'individual_variability', 'voxelwise_variability.mat'), 'voxelwise_variability');
load(fullfile('data', 'init.mat'), 'initV');

% Output figure directory
output_dir = fullfile('output', 'individual_variability');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

NetNum = 11;
[~, hard_labels] = max(initV, [], 2);  % Assign each voxel to max-loading network


%% ------------------- Group Variability by Network ---------------------
network_variability = cell(NetNum, 1);
network_index = cell(NetNum, 1);

for net_id = 1:NetNum
    vox_idx = find(hard_labels == net_id);
    network_variability{net_id} = voxelwise_variability(vox_idx)';
    network_index{net_id} = net_id * ones(1, length(vox_idx));
end

% Compute network-wise median to rank networks by variability
median_vals = zeros(NetNum, 1);
for net_id = 1:NetNum
    median_vals(net_id) = median(network_variability{net_id});
end
[~, sort_idx] = sort(median_vals);

% Re-assign x-axis positions based on rank
for rank = 1:NetNum
    network_index{sort_idx(rank)} = rank * ones(size(network_index{sort_idx(rank)}));
end

% Flatten for plotting
y = horzcat(network_variability{:});
x = horzcat(network_index{:});


%% ------------------- Define Colors ---------------------
% Primary vs. Association network color coding
color_primary = [162, 20, 47] / 255;
color_assoc   = [237, 177, 32] / 255;

C = [color_primary; color_primary; color_assoc; color_primary; ...
     color_assoc; color_primary; color_assoc; color_assoc; ...
     color_assoc; color_assoc; color_assoc];

%% ------------------- Plot Violin ---------------------
figureUnits = 'centimeters';
figureWidth = 15;
figureHeight = 6;

figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);
hold on;

violinplot(y, x, 'ViolinColor', {C}, 'ShowData', false, 'ViolinAlpha', 0.9);

ylabel('Across-subject variability', 'FontName', 'Arial', 'FontSize', 14);

set(gca, ...
    'Box', 'off', ...
    'LineWidth', 1, ...
    'TickDir', 'out', ...
    'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1], ...
    'FontName', 'Arial', 'FontSize', 14, ...
    'YTick', 0:0.02:1, ...
    'YLim', [0 0.09], ...
    'XTick', 1:NetNum, ...
    'XTickLabel', {'1','3','8','4','5','2','7','9','11','6','10'}, ...
    'YGrid', 'off');

set(gcf, 'Color', [1 1 1]);


%% ------------------- Export Figure ---------------------
set(figureHandle, 'PaperUnits', figureUnits);
set(figureHandle, 'PaperPosition', [0 0 figureWidth figureHeight]);

outname = fullfile(output_dir, 'individual_variability_network_violinplot.png');
print(figureHandle, outname, '-r300', '-dpng');
