%% visualize_group_networks.m 
% Visualize the 11 group-level soft parcellation networks using BrainNet Viewer
% Author: jlzhao
% Project: Personalized functional topography in neonatal brain networks

clc; clear;

%% ------------------------ Add BrainNet Viewer ------------------------
addpath(genpath('utils/BrainNetViewer'));

%% ------------------------ Path Settings ------------------------
output_dir        = fullfile('output', 'fig_grp_soft');
template_surface  = fullfile('surface_template', 'infant6_dHCP_V2.nv');
group_icn_dir     = fullfile('data', 'group_icn_soft');
color_map         = fullfile('colormap', 'green_yellow_red_deepred_custom.mat');

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% ------------------------ Network Names ------------------------
network_names = {'FroLim', 'aDMN', 'DorFro', 'FootMot', 'HandMot', ...
                 'SupPar', 'rDMN', 'lDMN', 'pDMN', 'LatVis', 'PriVis'};

%% ------------------------ Visualization Loop ------------------------
for i = 1:length(network_names)
    net_name = network_names{i};

    % Input NIfTI file
    if i < 10
        atlas_file = fullfile(group_icn_dir, ['icn_00', num2str(i), '.nii']);
    else
        atlas_file = fullfile(group_icn_dir, ['icn_0', num2str(i), '.nii']);
    end

    % Output image file named with network name
    fig_name = fullfile(output_dir, sprintf('GroupNetwork_%s.jpg', net_name));

    % Call BrainNet Viewer
    BrainNet_MapCfg(template_surface, atlas_file, color_map, fig_name);
end


