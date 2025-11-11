%% predict_brain_age.m
% Predict postmenstrual age using Ti-PCA + SVR framework
% Author: jlzhao

clear; clc; tic;

%% --------------------------- Configuration ---------------------------
addpath(genpath('utils/SVR'))  % Add SVR utilities to path

% Parameters
n_subjects = 315;
n_networks = 11;
n_voxels_per_net = 6905;
% network order: ['FroLim' 'aDMN' 'DorFro' 'FootMot' 'HandMot' 
% 'SupPar' 'rDMN' 'lDMN' 'pDMN' 'LatVis' 'PriVis']
network_order = [1,2,3,4,5,6,7,8,9,10,11];
reduction_method = 'Pca';
fold_quantity = 10;
pre_method = 'Normalize';
c_para = 1;  % SVR C value
perm_flag = 0;  % 0 for normal run, 1 for permutation test

%% --------------------------- Load Data ---------------------------
data_dir = 'data';
atlas_file = fullfile(data_dir, 'network_loading', 'V_vector_age_regressed.mat');
age_file   = fullfile(data_dir, 'PredictionAgeScore_Regressed.mat');

load(atlas_file);         % V_vector_age
load(age_file);           % PredictionAgeScore

SubjectsData = V_vector_age;

clear V_vector_age PredictionAgeScore;

%% --------------------------- Prediction ---------------------------
result_dir = fullfile('output', 'prediction', 'SVR_10CV_TiPCA_Age');
if ~exist(result_dir, 'dir'); mkdir(result_dir); end

Prediction = SVR_NFolds_Topography_PCA(...
    SubjectsData, Age, fold_quantity, pre_method, c_para, ...
    n_networks, n_voxels_per_net, ...
    reduction_method, network_order, perm_flag, result_dir);

save(fullfile(result_dir, 'Prediction.mat'), 'Prediction');

%% --------------------------- Evaluation ---------------------------
TrueScore = vertcat(Prediction.TrueScore{:});
PreScore = vertcat(Prediction.Score{:});
[r_value, p_value] = corr(TrueScore, PreScore);
save(fullfile(result_dir, 'Correlation.mat'), 'r_value', 'p_value');
disp(['Pearson correlation: r = ', num2str(r_value), ', p = ', num2str(p_value)]);
toc;


%% --------------------------- Compute Feature Weights ---------------------------
% This part calculates SVR contribution weights (not per fold)
weight_dir = fullfile(result_dir, 'weights');
if ~exist(weight_dir, 'dir'); mkdir(weight_dir); end

W_Calculate_SVR(SubjectsData, Age, pre_method, c_para, weight_dir);
