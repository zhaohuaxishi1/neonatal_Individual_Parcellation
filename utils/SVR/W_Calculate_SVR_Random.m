function [w_Brain, model_All] = W_Calculate_SVR_Random(Time, Subjects_Data, Subjects_Scores, FoldQuantity, Pre_Method, C_Range,Netnum, Netfea,randNet, ResultantFolder)
% Estimation of SVR contribution weights using repeated random
% stratified cross-validation.
%
% This function was used to identify high-contributing networks by
% quantifying the contribution weights of imaging features within
% repeated random 10-fold cross-validation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%




if ~exist(ResultantFolder, 'dir')
    mkdir(ResultantFolder);
end


[Subjects_Quantity, Features_Quantity] = size(Subjects_Data);

% Split into N folds randomly by bin
EachPart_Quantity = fix(Subjects_Quantity / FoldQuantity);
[~, SortedID] = sort(Subjects_Scores);
for j = 1:EachPart_Quantity
    BinID{j} =  SortedID((j - 1) * FoldQuantity + 1: j * FoldQuantity);
end

Origin_ID = cell(FoldQuantity, 1);
for j = 1:FoldQuantity
    for k = 1:EachPart_Quantity
        random_index = randi(length(BinID{k}));
        Origin_ID{j} = [Origin_ID{j} ;BinID{k}(random_index)];
        BinID{k}(random_index) = [];
    end
end


Reamin = mod(Subjects_Quantity, FoldQuantity);
ReaminBinID = SortedID(FoldQuantity * EachPart_Quantity + 1:end);
for j = 1:Reamin
    random_index = randi(length(ReaminBinID));
    Origin_ID{j} = [Origin_ID{j} ;ReaminBinID(random_index)];
    ReaminBinID(random_index) = [];
end

clear ReaminBinID random_index k BinID


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
    
    NetfeaPCA = Subjects_Quantity - EachPart_Quantity - 2;
    if mod(Subjects_Quantity, FoldQuantity) == 0
        NetfeaPCA = Subjects_Quantity - EachPart_Quantity - 1;
    end
    % data reduction
    Training_data_voxel = [];
    Test_data_voxel = [];
    for m = 1:Netnum
        % Traing test
        Training_data_voxel = [Training_data_voxel, Training_data(:,(randNet(m)-1)*Netfea+1:randNet(m)*Netfea)];
        Test_data_voxel = [Test_data_voxel,test_data(:,(randNet(m)-1)*Netfea+1:randNet(m)*Netfea)];
    end
      
    C_Optimal = C_Range;
    
    
    Training_data_New = Training_data_voxel;
    
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
    Training_data_final = double(Training_data_New);
    model_All = svmtrain(Training_scores, Training_data_final, ['-s 3 -t 0 -c ' num2str(C_Optimal)]);
    %     model = svmtrain(Training_scores, Training_data_final, ['-s 3 -t 2 -c ' num2str(C_Parameter)]);
    
    w_Brain = zeros(1, Features_Quantity);
    for k = 1 : model_All.totalSV
        w_Brain = w_Brain + model_All.sv_coef( k) * model_All.SVs(k, :);
    end

    w_Brain = w_Brain / norm(w_Brain);

    save([ResultantFolder filesep  'Prediction_' num2str(Time) '_w_Brain_fold_' num2str(fold) '.mat'], 'w_Brain');
end
