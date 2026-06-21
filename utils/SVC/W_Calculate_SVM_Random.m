
function [w_Brain, model_All] = W_Calculate_SVM_Random(Subjects_Data, Subjects_Label, Fold_Quantity, Pre_Method, ResultantFolder)

if nargin >= 5
    if ~exist(ResultantFolder, 'dir')
        mkdir(ResultantFolder);
    end
end

[Subjects_Quantity, Feature_Quantity] = size(Subjects_Data);
EachPart_Quantity = fix(Subjects_Quantity / Fold_Quantity);
[Splited_Data, Splited_Data_Label, Origin_ID_Cell] = Split_NFolds(Subjects_Data, Subjects_Label, Fold_Quantity);

predicted_labels = [];
decision_values = [];

for i = 1:Fold_Quantity
    
    disp(['The ' num2str(i) ' iteration!']);
    fold = i;
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
    
    %     W_Calculate_SVM(Training_all_data, Label, Pre_Method, ResultantFolder);
    Training_all_data = double(Training_all_data);
    model_All = svmtrain(Label, Training_all_data, '-t 0');
    w_Brain = zeros(1, Feature_Quantity);
    for j = 1 : model_All.totalSV
        w_Brain = w_Brain + model_All.sv_coef(j) * model_All.SVs(j, :);
    end
    w_Brain = w_Brain / norm(w_Brain);
    
    save([ResultantFolder filesep '_w_Brain_fold_' num2str(fold) '.mat'], 'w_Brain');
end