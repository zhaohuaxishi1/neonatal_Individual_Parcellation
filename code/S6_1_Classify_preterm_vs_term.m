%% classify_preterm_vs_term.m
% Classify preterm vs term-born neonates using Ti-PCA + SVM
% Author: jlzhao

clear; clc;

%% ------------------------ Configuration ------------------------
addpath(genpath('utils/SVC'));  % SVM utilities path

n_networks = 11;
n_voxels_per_net = 6905;
rand_network_order = 1:11;

result_dir = fullfile('output', 'classification', 'SVM_10CV_TiPCA');
if ~exist(result_dir, 'dir'); mkdir(result_dir); end

%% ------------------------ Load Data ------------------------
data_dir = 'data/classification';
feature_file = fullfile(data_dir, 'V_vector_FDAll.mat');
label_file = fullfile(data_dir, 'ClassificationlabelFDAll.mat');

load(feature_file);    % V_vector
load(label_file);      % Classificationlabel

SubjectsData = V_vector;
SubjectsLabel = Classificationlabel.label';

% Covariates: scan age, sex, mean FD
scan_age = Classificationlabel.scan_age;
sex = cellfun(@(x) double(x), Classificationlabel.sex(:,2));
meanFD = Classificationlabel.MeanFD;
covariates = [scan_age, sex, meanFD];

% Regress covariates from features
for i = 1:size(SubjectsData,2)
    [b, ~, res] = regress(SubjectsData(:,i), [ones(size(covariates,1),1), covariates]);
    SubjectsData(:,i) = res;
end

clear V_vector Classificationlabel b res

%% ------------------------ SVM Classification ------------------------
rng(89);  % Seed for reproducibility

fold_quantity = 10;
pre_method = 'Normalize';
perm_flag = 0;  % 0: normal run; 1: permutation test

[accuracy, sensitivity, specificity] = SVM_2group_NFolds_TiPCA(...
    SubjectsData, SubjectsLabel, ...
    rand_network_order, n_networks, n_voxels_per_net, ...
    fold_quantity, pre_method, result_dir, perm_flag);

% Save performance
save(fullfile(result_dir, 'Performance.mat'), 'accuracy', 'sensitivity', 'specificity');

%% ------------------------ Feature Contribution ------------------------
W_Calculate_SVM(SubjectsData, SubjectsLabel, pre_method, result_dir);
