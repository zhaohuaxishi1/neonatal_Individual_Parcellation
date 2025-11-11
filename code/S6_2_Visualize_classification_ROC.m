%% plot_preterm_vs_term_ROC.m
% Plot ROC curve and calculate AUC for Ti-PCA + SVM classification
% Author: jlzhao

clc; clear;

%% ------------------------ Configuration ------------------------
result_dir = fullfile('output', 'classification', 'SVM_10CV_TiPCA');
label_file = fullfile('data', 'classification', 'ClassificationlabelFDAll.mat');

addpath(genpath('utils/SVC'));  % ROC utility assumed to be here

%% ------------------------ Load Data ------------------------
% Decision values saved during classification
load(fullfile(result_dir, 'Y.mat'));  % Contains Y_group1, Y_group0
DecisionValues = [Y_group1, Y_group0]';  % N x 1

% Ground truth labels: 1 = term, -1 = preterm
load(label_file);  % Contains Classificationlabel.label
Label = Classificationlabel.label;  % N x 1

%% ------------------------ Plot ROC & Compute AUC ------------------------
ROC_Draw_Flag = 1;  % Set to 1 to display plot
AUC_Calculate_ROC_Draw(DecisionValues, Label, ROC_Draw_Flag,result_dir);

