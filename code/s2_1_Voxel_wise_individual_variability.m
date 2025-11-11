%% compute_visualize_individual_variability.m
% Compute and visualize individual variability based on probabilistic network maps
% Author: jlzhao

clc; clear;

%% Step 1: Compute variability map using MAD
% Load individual probabilistic network decompositions
addpath(genpath('utils'));

load('data/individual_parcellation/IndividualAtlas.mat');  % contains final_UV
 
num_networks = 11;
V_matrix = [];

% Concatenate all individual V matrices
for subj = 1:size(final_UV, 1)
    disp(['Processing subject ', num2str(subj)]);
    individual_V = final_UV{subj, 1};
    V_matrix = [V_matrix, individual_V{1}];
end

% Compute MAD across subjects per voxel for each network
V_network = [];
for k = 1:num_networks
    net_k_matrix = V_matrix(:, k:num_networks:end);  % extract current network slice
    V_network = [V_network; mad(net_k_matrix', 1)];
end

% Mean MAD across networks for each voxel
voxelwise_variability = mean(V_network)';

% Save results
mkdir('output/individual_variability');
save('output/individual_variability/voxelwise_variability.mat', 'voxelwise_variability');

% Save raw matrix for visualization or reordering
V_network = V_network';
save('output/individual_variability/V_network.mat', 'V_network');


%% Step 2: Save variability map as NIfTI for visualization
% Prepare image projection
% Convert to NIfTI
refNii = 'data/Resliced3mm_GmMask_final.nii';

func_saveVol2Nii('voxelwise_variability', ...
    'output/individual_variability/voxelwise_variability.mat', ...
    refNii, ...
    'output/individual_variability', ...
    'Ind_Var_Map');

%% Step 3: Visualization using BrainNet Viewer

% Add BrainNet Viewer path
addpath(genpath('utils/BrainNetViewer'));

% Set visualization paths
output_dir       = 'output/individual_variability';
nii_path         = fullfile(output_dir, 'Ind_Var_Map.nii');  % this is the saved NIfTI file
surface_template = 'surface_template/infant6_dHCP_V2.nv';
color_map        = 'colormap/IndVar_colormap.mat';
fig_path         = fullfile(output_dir, 'Ind_Var_BrainMap.jpg');

% Generate brain map
BrainNet_MapCfg(surface_template, nii_path, color_map, fig_path);



