function Prediction = SVR_NFolds_Topography_PCA_RobustCovariate(Subjects_Data, Subjects_Scores,cov_new, FoldQuantity, Pre_Method, C_Range,Netnum, Netfea, Reduction_Method,randNet,Permutation_Flag, ResultantFolder)
%% ------------------------------------------------------------------------
% SVR Prediction with Robust Covariate Regression
%
% Author: Jianlong Zhao
%
% Description
% ------------------------------------------------------------------------
% Sensitivity analysis using robust covariate regression.
%
% All participants are retained in the prediction framework. Covariate
% effects are removed using robust regression (robustfit) instead of
% ordinary least-squares regression.
%
% To avoid information leakage, covariate regression is performed
% exclusively within each training fold and subsequently applied to
% the corresponding test fold.

if ~exist(ResultantFolder,'dir')
    mkdir(ResultantFolder);
end

[Subjects_Quantity, ~] = size(Subjects_Data);

% Split into N folds randomly
EachPart_Quantity = fix(Subjects_Quantity / FoldQuantity);

%% Example:
ResPath = '...';

load(fullfile(ResPath,'Prediction.mat'));
load(fullfile(ResPath,'idx.mat'));

Subjects_Scores_del = Subjects_Scores(idx);

[~,SortedID_del] = sort(Subjects_Scores_del);

idx_rank = idx(SortedID_del)';

Origin_ID_top = Prediction.Origin_ID_top;
Origin_ID = Prediction.Origin_ID;

clear Prediction

specialFold = [7 10];

for f = 1:2
    Origin_ID_top{specialFold(f)} = ...
        [Origin_ID_top{specialFold(f)};
        idx_rank(f:2)];
end

%% ========================= Cross Validation ===========================
for j = 1:FoldQuantity
    
    disp(['The ' num2str(j) ' fold!']);
    
    if  ~ismember(j,index)
        mask = true(size(Subjects_Data,1),1);
        mask(idx) = false;
        Training_data = Subjects_Data(mask,:);
        Training_scores = Subjects_Scores(mask,:);
        
        test_data = Training_data(Origin_ID{j}, :);
        test_score = Training_scores(Origin_ID{j});
        Training_data(Origin_ID{j}, :) = [];
        Training_scores(Origin_ID{j}) = [];
        Training_cov_new = cov_new(mask,:);
        test_cov_new = Training_cov_new(Origin_ID{j},:);
        Training_cov_new(Origin_ID{j},:) = [];
        
    else
        Training_data = Subjects_Data;
        Training_scores = Subjects_Scores;
        % Select training data and testing data
        test_data = Training_data(Origin_ID_top{j}, :);
        test_score = Training_scores(Origin_ID_top{j});
        Training_data(Origin_ID_top{j}, :) = [];
        Training_scores(Origin_ID_top{j}) = [];
        
        Training_cov_new = cov_new;
        Training_cov_new(Origin_ID_top{j},:) = [];
        test_cov_new = cov_new(Origin_ID_top{j},:);
    end
    
    % demean
    covariates_mean = mean(Training_cov_new,1);
    Training_cov_new = Training_cov_new-repmat(covariates_mean,[length(Training_scores) 1]);
    test_cov_new = test_cov_new-repmat(covariates_mean,[length(test_score) 1]);
    %% Regressing covariates from behavior
    [yb,~] = robustfit(Training_cov_new, Training_scores);
    Training_data_score = Training_scores-Training_cov_new*yb(2:end);%
    Test_data_score = test_score-test_cov_new*yb(2:end);%
    
    Training_scores = Training_data_score;
    test_score = Test_data_score;
    clear Training_data_score Test_data_score
    
    
    %%
    if Permutation_Flag
        Training_Quantity = length(Training_scores);
        RandIndex = randperm(Training_Quantity);
        Training_scores = Training_scores(RandIndex);
    end
    
    %% Normalize Training data
    Training_data_New = Training_data;
    
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
    
    Training_data  = Training_data_New;
    
    %% Normalize test data
    test_data_New = test_data;
    
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
    
    test_data  = test_data_New;
    %%
    
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
    
    clear Indicator_Tmp Training_data_voxel_pca Training_data_voxel Test_data_voxel_pca Test_data_voxel;
    
    for i = 1:size(Training_data_pca,2)
        [b(:,i),~] = robustfit(Training_cov_new, Training_data_pca(:,i));
        Training_data_pca(:,i) = Training_data_pca(:,i)-Training_cov_new*b(2:end,i);
        Test_data_pca(:,i) = Test_data_pca(:,i)-test_cov_new*b(2:end,i);
    end
    
    %%
    C_Optimal = C_Range;
    
    
    %%  Training Normalize
    Training_data_New = Training_data_pca;
    %%
    % SVR training
    Training_data_final = double(Training_data_New);
    model = svmtrain(Training_scores, Training_data_final, ['-s 3 -t 0 -c ' num2str(C_Optimal)]);
    
    %%
    test_data_New = Test_data_pca;
    clear Test_data_pca Training_data_pca
    test_data_final = double(test_data_New);
    
    
    %% Predict test data
    [Predicted_Scores, ~, ~] = svmpredict(test_score, test_data_final, model);
    Prediction.Origin_ID{j} = Origin_ID{j};
    Prediction.TrueScore{j} = test_score;
    Prediction.Score{j} = Predicted_Scores;
    Prediction.Corr(j) = corr(Predicted_Scores, test_score);
    Prediction.y_regcoff(:,j) = yb;
    
    corr(Predicted_Scores, test_score)
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

