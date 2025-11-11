%% visualize_group_networks.m
% Visualize the 11 group-level functional networks using BrainNet Viewer
% Author: jlzhao
clc; clear;
%% Add BrainNet Viewer (if needed)
addpath(genpath('utils/BrainNetViewer'));
%% Initialization
% Set paths
output_dir = fullfile('output', 'fig_grp_hard');
template_surface = fullfile('surface_template', 'infant6_dHCP_V2.nv');
group_icn_dir = fullfile('data', 'group_icn_hard');
color_map = fullfile('colormap', 'Hard_atlas.mat');

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


%% Visualization Loop

atlas_file = fullfile(group_icn_dir, 'group_hard.nii');

% Output figure name
fig_name =[output_dir '\HardAtlas.jpg'];

% Call BrainNet Viewer
BrainNet_MapCfg(template_surface, atlas_file, color_map, fig_name);


