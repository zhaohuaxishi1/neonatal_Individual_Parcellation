function [Accuracy, Sensitivity, Specificity] = SVM_2group_NFolds_TiPCA(Subjects_Data, Subjects_Label,randNet,Netnum,Netfea, Fold_Quantity, Pre_Method,ResultantFolder, Permutation_Flag)
% SVM_2group_NFolds_TiPCA
% -------------------------------------------------------------------------
% Perform N-fold SVM classification with Ti-PCA (per-network PCA) 
% on multi-network feature matrices. Typically used to classify two groups
% (e.g., preterm vs term infants) using topographic features.
%
% ------------------------- Input Arguments -------------------------------
% Subjects_Data     - [m x n] matrix
%                    m = number of subjects, n = number of features (across networks)
%
% Subjects_Label    - [m x 1] vector
%                    Binary class labels (-1 or 1) for each subject
%
% randNet           - [1 x Netnum] vector
%                    Order of networks to be processed (index-based)
%
% Netnum            - Integer, number of functional networks (e.g., 11)
%
% Netfea            - Integer, number of voxels/features per network
%
% Fold_Quantity     - Integer, number of folds for cross-validation
%
% Pre_Method        - String, preprocessing method: 'Normalize' or 'Scale'
%
% ResultantFolder   - Output directory to store classification results
%
% Permutation_Flag  - 0: regular CV; 1: permutation test (shuffles labels)
%
% ------------------------- Output Arguments ------------------------------
% Accuracy          - Overall classification accuracy (float)
% Sensitivity       - True positive rate (class -1)
% Specificity       - True negative rate (class 1)
%
% ------------------------- Outputs Saved (if Permutation_Flag is set) ----
%   - Y.mat:           decision values for each group
%   - Category.mat:    predicted labels for each group
%   - Accuracy.mat:    accuracy value
%   - Sensitivity.mat: sensitivity value
%   - Specificity.mat: specificity value
%   - WrongInfo.mat:   misclassified subject IDs and counts
%   - Feature weights are computed and saved via `W_Calculate_SVM`
%
% Author: jlzhao
% -------------------------------------------------------------------------


if nargin > 8 & ~exist(ResultantFolder, 'dir')
    mkdir(ResultantFolder);
end

[Subjects_Quantity, Feature_Quantity] = size(Subjects_Data);
EachPart_Quantity = fix(Subjects_Quantity / Fold_Quantity);
[Splited_Data, Splited_Data_Label, Origin_ID_Cell] = Split_NFolds(Subjects_Data, Subjects_Label, Fold_Quantity);

predicted_labels = [];
decision_values = [];
for i = 1:Fold_Quantity
    
    disp(['The ' num2str(i) ' iteration!']);
    
    % Select training data and testing data
    test_label = Splited_Data_Label{i};
    test_data = Splited_Data{i};
    
    Training_all_data = [];
    Label = [];
    for j = 1:Fold_Quantity
        if j == i
            continue;
        end
        Training_all_data = [Training_all_data; Splited_Data{j}];
        Label = [Label; Splited_Data_Label{j}];
    end
    
    if Permutation_Flag
        Rand_ID = randperm(length(Label));
        Label = Label(Rand_ID);
    end
    
    if strcmp(Pre_Method, 'Normalize')
        %Normalizing
        MeanValue = mean(Training_all_data);
        StandardDeviation = std(Training_all_data);
        [rows, columns_quantity] = size(Training_all_data);
        for j = 1:columns_quantity
            Training_all_data(:, j) = (Training_all_data(:, j) - MeanValue(j)) / StandardDeviation(j);
        end
    elseif strcmp(Pre_Method, 'Scale')
        % Scaling to [0 1]
        MinValue = min(Training_all_data);
        MaxValue = max(Training_all_data);
        [rows, columns_quantity] = size(Training_all_data);
        for j = 1:columns_quantity
            Training_all_data(:, j) = (Training_all_data(:, j) - MinValue(j)) / (MaxValue(j) - MinValue(j));
        end
    end
    Training_all_data(isnan(Training_all_data)==1) = 0;
    
    if strcmp(Pre_Method, 'Normalize')
        % Normalizing
        MeanValue_New = repmat(MeanValue, length(test_label), 1);
        StandardDeviation_New = repmat(StandardDeviation, length(test_label), 1);
        test_data = (test_data - MeanValue_New) ./ StandardDeviation_New;
    elseif strcmp(Pre_Method, 'Scale')
        % Scale
        MaxValue_New = repmat(MaxValue, length(test_label), 1);
        MinValue_New = repmat(MinValue, length(test_label), 1);
        test_data = (test_data - MinValue_New) ./ (MaxValue_New - MinValue_New);
    end
    test_data(isnan(test_data)==1) = 0;
    
    
    NetfeaPCA = Subjects_Quantity - EachPart_Quantity - 3;
    if mod(Subjects_Quantity, Fold_Quantity) == 0
        NetfeaPCA = Subjects_Quantity - EachPart_Quantity - 2;
    end
    % data reduction
    Training_data_pca = [];
    Test_data_pca = [];
    coeff_pca =  [];
    for m = 1:Netnum
        % Traing test
        Training_data_voxel = Training_all_data(:,(randNet(m)-1)*Netfea+1:randNet(m)*Netfea);
        Test_data_voxel = test_data(:,(randNet(m)-1)*Netfea+1:randNet(m)*Netfea);
        
        [coeff,score,latent,tsquared,explained,~] = pca(Training_data_voxel,'NumComponents',NetfeaPCA);
        Training_data_voxel_pca = Training_data_voxel * coeff(:,:);
        Training_data_pca = [Training_data_pca,Training_data_voxel_pca];
        coeff_pca =  [coeff_pca;coeff];
        % Testing test
        Test_data_voxel_pca = Test_data_voxel * coeff(:,:);
        Test_data_pca = [Test_data_pca,Test_data_voxel_pca];
        
    end
    
    clear Indicator_Tmp Training_data_voxel_pca Training_data_voxel Test_data_voxel_pca Test_data_voxel;
    
    Training_all_data = Training_data_pca;
    
    
    % classification
    Label = reshape(Label, length(Label), 1);
    Training_all_data = double(Training_all_data);
    model(i) = svmtrain(Label, Training_all_data, '-t 0');
    
    
    test_data = Test_data_pca;
    % predicts
    test_data = double(test_data);
    [predicted_labels_tmp, ~, ~] = svmpredict(test_label, test_data, model(i));
    predicted_labels = [predicted_labels predicted_labels_tmp'];
        
    ComponentQuantity = size(test_data,2);
    w{i} = zeros(1, ComponentQuantity);
    for j = 1 : model(i).totalSV
        w{i} = w{i} + model(i).sv_coef(j) * model(i).SVs(j, :);
    end
    decision_values_tmp = w{i} * test_data' - model(i).rho;
    decision_values = [decision_values decision_values_tmp];
    
end

Origin_ID = [];
for i = 1:length(Origin_ID_Cell)
    Origin_ID = [Origin_ID; Origin_ID_Cell{i}];
end

Group1_Index = find(Subjects_Label(Origin_ID) == 1);
Group0_Index = find(Subjects_Label(Origin_ID) == -1);
Category_group1 = predicted_labels(Group1_Index);
Y_group1 = decision_values(Group1_Index);
Category_group0 = predicted_labels(Group0_Index);
Y_group0 = decision_values(Group0_Index);

group0_Wrong_ID = find(Category_group0 == 1);
group0_Wrong_Quantity = length(group0_Wrong_ID);
group1_Wrong_ID = find(Category_group1 == -1);
group1_Wrong_Quantity = length(group1_Wrong_ID);
disp(['group0: ' num2str(group0_Wrong_Quantity) ' subjects are wrong ' mat2str(group0_Wrong_ID) ]);
disp(['group1: ' num2str(group1_Wrong_Quantity) ' subjects are wrong ' mat2str(group1_Wrong_ID) ]);
Accuracy = (Subjects_Quantity - group0_Wrong_Quantity - group1_Wrong_Quantity) / Subjects_Quantity;
disp(['Accuracy is ' num2str(Accuracy) ' !']);
Sensitivity = (length(Group0_Index) - group0_Wrong_Quantity) / length(Group0_Index);
Specificity = (length(Group1_Index) - group1_Wrong_Quantity) / length(Group1_Index);
disp(['Sensitivity is ' num2str(Sensitivity) ' !']);
disp(['Specificity is ' num2str(Specificity) ' !']);

if nargin > 8
    save([ResultantFolder filesep 'Y.mat'], 'Y_group1', 'Y_group0');
    save([ResultantFolder filesep 'Category.mat'], 'Category_group1', 'Category_group0');
    save([ResultantFolder filesep 'WrongInfo.mat'], 'group0_Wrong_Quantity', 'group0_Wrong_ID', 'group1_Wrong_Quantity', 'group1_Wrong_ID');
    save([ResultantFolder filesep 'Accuracy.mat'], 'Accuracy');
    save([ResultantFolder filesep 'Sensitivity.mat'], 'Sensitivity');
    save([ResultantFolder filesep 'Specificity.mat'], 'Specificity');
    % Calculating w
    W_Calculate_SVM(Subjects_Data, Subjects_Label, Pre_Method, ResultantFolder);
end

