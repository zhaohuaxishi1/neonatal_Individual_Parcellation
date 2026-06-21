%% predict_language_score.m
% Predict 18-month language outcome using Ti-PCA + SVR framework
% Author: jlzhao

clear; clc; tic;

%% --------------------------- Configuration ---------------------------
addpath(genpath('utils/SVR'));  % Add SVR utilities to path

% Parameters
n_subjects = 249;
n_networks = 11;
n_voxels_per_net = 6905;
% network order: ['FroLim' 'aDMN' 'DorFro' 'FootMot' 'HandMot' 
% 'SupPar' 'rDMN' 'lDMN' 'pDMN' 'LatVis' 'PriVis']
network_order = [1,2,3,4,5,6,7,8,9,10,11];
reduction_method = 'Pca';
fold_quantity = 10;
pre_method = 'Normalize';
c_para = 1;
perm_flag = 0;

% Paths
data_dir = 'data';
atlas_file = fullfile(data_dir, 'network_loading', 'V_vector_cog_prediction.mat');
score_file = fullfile(data_dir, 'PredictionLanguageScore.mat');

result_dir = fullfile('output', 'prediction', 'SVR_10CV_TiPCA_Language');
if ~exist(result_dir, 'dir'); mkdir(result_dir); end

%% --------------------------- Load Data ---------------------------
load(atlas_file);         % -> V_vector_age
SubjectsData = V_vector;

load(score_file);         % -> PredictionCognitionScore

clear V_vector PredictionCognitionScore;

%% --------------------------- Prediction ---------------------------
Prediction = SVR_NFolds_Topography_PCA(...
    SubjectsData, Language, fold_quantity, pre_method, c_para, ...
    n_networks, n_voxels_per_net, ...
    reduction_method, network_order, perm_flag, result_dir);

save(fullfile(result_dir, 'Prediction.mat'), 'Prediction');

%% --------------------------- Evaluation ---------------------------
TrueScore = vertcat(Prediction.TrueScore{:});
PreScore  = vertcat(Prediction.Score{:});
[r_value, p_value] = corr(TrueScore, PreScore);
MAE = mean(abs(TrueScore - PreScore));
save(fullfile(result_dir, 'Correlation.mat'), 'r_value', 'p_value');

fprintf('[Language Prediction] r = %.3f, p = %.3g\n', r_value, p_value);
toc;
%% --------------------------- Compute Feature Weights ---------------------------
weight_dir = fullfile(result_dir, 'weights');
if ~exist(weight_dir, 'dir'); mkdir(weight_dir); end

for i = 1: 100
    W_Calculate_SVR_Random(i, SubjectsData, Language, fold_quantity, pre_method, c_para, n_networks, n_voxels_per_net, network_order, weight_dir);
end

