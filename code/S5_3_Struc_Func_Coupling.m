% PLS Analysis of Structural Maturation and Functional Topographic Refinement
%
% Author: Jianlong Zhao
%
% This script evaluates the relationship between cortical structural
% maturation and functional topographic refinement using Partial Least
% Squares (PLS) analysis.
%
% To mitigate potential biases arising from unequal network sizes,
% a repeated size-matched subsampling procedure was performed prior
% to PLS analysis. Specifically, 190 vertices were randomly selected
% from each network, and the analysis was repeated 100 times to assess
% the robustness of latent variable correlations and feature loadings.
%
% Analysis:
%   1. Whole-brain PLS analysis
%   2. Network-level PLS analysis
%
% Inputs:
%   AgeStrufeature.mat
%   Agefuncfeature.mat
%   group_data.mat
%
% Outputs:
%   Bootstrap PLS results
%% ------------------------------------------------------------------------
clc;
clear;
addpath('utils/PLS_Codes');
addpath('utils');

%% Paths
data_dir = 'data';
% --------------------------- Load Features ---------------------------
load(fullfile(data_dir,'AgeStrufeature.mat'));   % structural maturation features
load(fullfile(data_dir,'Agefuncfeature.mat'));   % functional topographic refinement features

result_dir = fullfile('output', 'PLS_Struc_Func_Age');
if ~exist(result_dir, 'dir'); mkdir(result_dir); end

%% --------------------------- Whole-Brain PLS ---------------------------
rng('default')
r_p = zeros(100,2);
overall_LC_pvals = zeros(100,1);
Str_contri = zeros(100,4);
numSamples = 190;

for n_boot = 1:100
    fprintf('Whole-brain bootstrap %d/%d\n',n_boot);
    [sampled_sturcture_feature,sampled_varfeatures,index] = random_uniform_sample(sturcture_feature, varfeature,numSamples);
    
    % Remove vertices with missing or zero values
    nonzero_index_A = find(sampled_sturcture_feature(:,1) ~= 0);
    nonzero_index_B = find(sampled_varfeatures(:,1) ~= 0);
    common_index = intersect(nonzero_index_A, nonzero_index_B);
    sampled_varfeatures = sampled_varfeatures(common_index,:);
    sampled_sturcture_feature = sampled_sturcture_feature(common_index,:);
    sampled_sturcture_feature(isnan(sampled_sturcture_feature)) = 0;
    
    X0 = sampled_sturcture_feature;
    
    % Run PLS analysis
    [input,pls_opts,save_opts] = myPLS_inputs(X0,sampled_varfeatures,10000,1000);
    [input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);
    res = myPLS_analysis(input,pls_opts);
    mkdir([result_dir '/Bootstrap']);
    save([result_dir '/Bootstrap/PLS_age_correlation_all_boots_' num2str(n_boot)],'res');
    
    [r_p(n_boot,1),r_p(n_boot,2)] = corr(res.Lx,res.Ly);
    overall_LC_pvals(n_boot,1) = res.LC_pvals;
    Str_contri(n_boot,:) = res.LC_img_loadings;
end

%% ------------------------------------------------------------------------
% Network-Level PLS Analysis
%
% Networks differ substantially in spatial extent. To ensure that
% PLS results are not driven by network size differences, networks
% containing more than 190 vertices were randomly downsampled to
% 190 vertices prior to each analysis.
%
% The procedure was repeated 100 times for each network.
%% ------------------------------------------------------------------------
clc;
clear;
addpath('D:\software\Software Package\PLS_Codes');

load('group_data');

nNetworks = 11;

network_length = zeros(nNetworks,1);

for i = 1:nNetworks
    network_length(i) = sum(group_data == i);
end

% --------------------------- Load Features ---------------------------
load(fullfile(data_dir,'AgeStrufeature.mat'));   % structural maturation features
load(fullfile(data_dir,'Agefuncfeature.mat'));   % functional topographic refinement features
result_dir = fullfile('output', 'PLS_Struc_Func_Age');
if ~exist(result_dir, 'dir'); mkdir(result_dir); end

%%
rng('default');
for n_boot = 1:100
    for i = 1:11
        fprintf('Bootstrap %d/%d | Network %d/%d\n', n_boot,i,nNetworks);
        index = find(group_data == i);
        indvar_network = varfeature(index,:);
        sturcture_network = sturcture_feature(index,:);
        numSamples = 190;
        
        % Downsample large networks to generate size-matched datasets
        if network_length(i,1) > numSamples
            [sampled_sturcture_feature,sampled_varfeatures,index] = random_uniform_sample(sturcture_network ,indvar_network,numSamples);
            X0=sampled_sturcture_feature;
            [input,pls_opts,save_opts] = myPLS_inputs(X0,sampled_varfeatures,10000,1000);
        else
            X0=sturcture_network;
            [input,pls_opts,save_opts] = myPLS_inputs(X0,indvar_network,10000,1000);
        end
        
        [input,pls_opts,save_opts] = myPLS_initialize(input,pls_opts,save_opts);
        res = myPLS_analysis(input,pls_opts);
        Str_contri(n_boot,:) = res.LC_img_loadings;
        
        [r_p(n_boot,i),~] = corr(res.Lx,res.Ly);
        overall_LC_pvals(n_boot,i) = res.LC_pvals;
        
        mkdir([result_dir '/Bootstrap']);
        save([result_dir '/Bootstrap/PLS_age_correlation_network_' num2str(i) '_boots_' num2str(n_boot)],'res')        
    end
end

