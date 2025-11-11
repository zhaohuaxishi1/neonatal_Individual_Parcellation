function Prediction = SVR_NFolds_Topography_PCA(Subjects_Data, Subjects_Scores, FoldQuantity, Pre_Method, C_Para,Netnum, Netfea,Reduction_Method,randNet,Permutation_Flag, ResultantFolder)
%% SVR_NFolds_Topography_PCA
% -------------------------------------------------------------------------
% Multivariate prediction framework combining topography-based feature modeling
% and within-network PCA dimensionality reduction before SVR regression.
%
% This implementation extends the pattern regression pipeline of 
% Zaixu Cui (Cui et al., 2016, 2018) by introducing topography-informed PCA 
% (Ti-PCA) within each functional network to enhance interpretability of
% network-specific representations.
%
% Author: jlzhao (based on original SVR_NFolds_Nest_PCA_CArray_Voxel_11Net by Zaixu Cui)
% Contact: jlzhao (Brain Function & Development Lab)
%
% References:
%   Cui et al., Cerebral Cortex, 2018
%   Cui & Gong, NeuroImage, 2018
%   Cui et al., Human Brain Mapping, 2016
% -------------------------------------------------------------------------
%
% INPUTS
%   Subjects_Data:     [m × n] matrix
%                      Each row is a subject, each column is a feature (voxel loading)
%
%   Subjects_Scores:   [m × 1] vector of target variable (e.g., postmenstrual age)
%
%   FoldQuantity:      Number of cross-validation folds (recommended: 5 or 10)
%
%   Pre_Method:        Preprocessing method for input features:
%                      'Normalize' | 'Scale' | 'None'
%
%   C_Para:            SVR regularization parameter (default = 1)
%
%   Netnum:            Number of functional networks
%
%   Netfea:            Number of voxels (features) per network
%
%   Reduction_Method:  Feature reduction method per network ('PCA')
%
%   randNet:           Index order of functional networks
%
%   Permutation_Flag:  1 = Perform permutation testing (shuffle target values)
%                      0 = Normal prediction
%
%   ResultantFolder:   Output folder for saving results (.mat)
%
% OUTPUT
%   Prediction:        Structure with fields:
%                      .TrueScore    — true labels per fold
%                      .Score        — predicted scores per fold
%                      .Corr         — Pearson r per fold
%                      .r_value      — overall correlation
%                      .p_value      — overall p-value
%                      .MAE          — mean absolute error
%
% -------------------------------------------------------------------------
% Example:
%   Prediction = SVR_NFolds_Topography_PCA(Data, Age, 10, 'Normalize', 1, ...
%                11, 6905, 'PCA', 1:11, 0, 'output/SVR_TiPCA');
% -------------------------------------------------------------------------
%
% Original framework:
%   Written by Zaixu Cui (zaixucui@gmail.com; Zaixu.Cui@pennmedicine.upenn.edu)
%   Modified and extended by jlzhao for individualized topographic PCA modeling
% -------------------------------------------------------------------------

if nargin >= 14
    if ~exist(ResultantFolder, 'dir')
        mkdir(ResultantFolder);
    end
end

[Subjects_Quantity, ~] = size(Subjects_Data);

% Split into N folds randomly
EachPart_Quantity = fix(Subjects_Quantity / FoldQuantity);
[~, SortedID] = sort(Subjects_Scores);
for j = 1:FoldQuantity
    Origin_ID{j} = SortedID([j : FoldQuantity : Subjects_Quantity]);
end

for j = 1:FoldQuantity
    
    disp(['The ' num2str(j) ' fold!']);
    fold = j;
    Training_data = Subjects_Data;
    Training_scores = Subjects_Scores;
    
    % Select training data and testing data
    test_data = Training_data(Origin_ID{j}, :);
    test_score = Training_scores(Origin_ID{j});
    Training_data(Origin_ID{j}, :) = [];
    Training_scores(Origin_ID{j}) = [];
    
    if Permutation_Flag
        Training_Quantity = length(Training_scores);
        RandIndex = randperm(Training_Quantity);
        Training_scores = Training_scores(RandIndex);
    end
    
    
    NetfeaPCA = Subjects_Quantity - EachPart_Quantity - 2;
    if mod(Subjects_Quantity, FoldQuantity) == 0
        NetfeaPCA = Subjects_Quantity - EachPart_Quantity - 1;
    end
    % data reduction
    Training_data_pca = [];
    Test_data_pca = [];
    coeff_pca =  [];
    for m = 1:Netnum
        % Traing test
        Training_data_voxel = Training_data(:,(randNet(m)-1)*Netfea+1:randNet(m)*Netfea);
        Test_data_voxel = test_data(:,(randNet(m)-1)*Netfea+1:randNet(m)*Netfea);
        if strcmp(Reduction_Method, 'Pca')
            [coeff,score,latent,tsquared,explained,~] = pca(Training_data_voxel,'NumComponents',NetfeaPCA);
            Training_data_voxel_pca = Training_data_voxel * coeff(:,:);
            Training_data_pca = [Training_data_pca,Training_data_voxel_pca];
            coeff_pca =  [coeff_pca;coeff];
            % Testing test
            Test_data_voxel_pca = Test_data_voxel * coeff(:,:);
            Test_data_pca = [Test_data_pca,Test_data_voxel_pca];
        end
    end
    %     NetfeaPCA = size(Training_data_voxel,1) + Indicatorfea - 1;
    
    clear Indicator_Tmp Training_data_voxel_pca Training_data_voxel Test_data_voxel_pca Test_data_voxel;
        
    C_Optimal = C_Para;
    
    
    Training_data_New = Training_data_pca;
    
    if strcmp(Pre_Method, 'Normalize')
        % Normalizing
        MeanValue = mean(Training_data_New);
        StandardDeviation = std(Training_data_New);
        [~, columns_quantity] = size(Training_data_New);
        for k = 1:columns_quantity
            if Training_data_New(:, k) == 0
                continue;
            else
                Training_data_New(:, k) = (Training_data_New(:, k) - MeanValue(k)) / StandardDeviation(k);
            end
        end
    elseif strcmp(Pre_Method, 'Scale')
        % Scaling to [0 1]
        MinValue = min(Training_data_New);
        MaxValue = max(Training_data_New);
        [~, columns_quantity] = size(Training_data_New);
        for k = 1:columns_quantity
            if Training_data_New(:, k) == 0
                continue;
            else
                Training_data_New(:, k) = (Training_data_New(:, k) - MinValue(k)) / (MaxValue(k) - MinValue(k));
            end
        end
    end
    
    % SVR training
%     Training_scores = Training_scores;
    Training_data_final = double(Training_data_New);
    model = svmtrain(Training_scores, Training_data_final, ['-s 3 -t 0 -c ' num2str(C_Optimal)]);
    
    test_data_New = Test_data_pca;
    % Normalize test data
    if strcmp(Pre_Method, 'Normalize')
        % Normalizing
        MeanValue_New = repmat(MeanValue, length(test_score), 1);
        StandardDeviation_New = repmat(StandardDeviation, length(test_score), 1);
        test_data_New = (test_data_New - MeanValue_New) ./ StandardDeviation_New;
        test_data_New(isnan(test_data_New)==1) = 0;
    elseif strcmp(Pre_Method, 'Scale')
        % Scale
        MaxValue_New = repmat(MaxValue, length(test_score), 1);
        MinValue_New = repmat(MinValue, length(test_score), 1);
        test_data_New = (test_data_New - MinValue_New) ./ (MaxValue_New - MinValue_New);
        test_data_New(isnan(test_data_New)==1) = 0;
    end
    test_data_final = double(test_data_New);
    % Predict test data
    [Predicted_Scores, ~, ~] = svmpredict(test_score, test_data_final, model);
    Prediction.Origin_ID{j} = Origin_ID{j};
    Prediction.TrueScore{j} = test_score;
    Prediction.Score{j} = Predicted_Scores;
    Prediction.Corr(j) = corr(Predicted_Scores, test_score);    
end

TrueScore = [];
PreScore = [];
for i = 1:FoldQuantity
    TrueScore = [TrueScore;Prediction.TrueScore{i}];
    PreScore = [PreScore;Prediction.Score{i}];
end

[r_value,p_value] = corr(TrueScore,PreScore);
Prediction.r_value = r_value;
Prediction.p_value = p_value;
Prediction.MAE = mean(abs(TrueScore - PreScore));
save([ResultantFolder filesep 'Prediction.mat'], 'Prediction');
