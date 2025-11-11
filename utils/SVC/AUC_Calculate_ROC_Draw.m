
function [AUC LargestAccuracy] = AUC_Calculate_ROC_Draw(DecisionValues, Label, ROC_Draw_Flag,result_dir)

%
% Decision_Values: n * 1
%                  y = wx + b
%                  if y > 0, positive class
%                  if y < 0, negative class
% 
% Label: 
%                  n * 1, each element is 1 or -1
%                  -1: patient
%                  1: NC
%
% ROC_Draw_Flag:
%                  1 or 0
%                  1: draw ROC curve
%                  0: don't draw ROC curve
%

%
% Written by Zaixu Cui, State Key Laboratory of Cognitive 
% Neuroscience and Learning, Beijing Normal University, 2013.
% Maintainer: zaixucui@gmail.com
%

[DecisionValues_rows, DecisionValues_columns] = size(DecisionValues);
if DecisionValues_columns ~= 1
    error('DecisionValues should be a n*1 vector!');
end
[Label_rows, Label_columns] = size(Label);
if Label_columns ~= 1
    error('Label should be a n*1 vector!');
end

P = length(find(Label == -1));
N = length(find(Label == 1));
TP = 0;
FP = 0;
TP_prev = 0;
FP_prev = 0;
[Sorted_DecisionValues OriginPos] = sort(DecisionValues, 1, 'ascend');
SubjectQuantity = length(Sorted_DecisionValues);

DecisionValue_prev = -1000000;
AUC = 0;

TP_Array = 0;
FP_Array = 0;
Accuracy_Array = N / (P + N);
for i = 1:SubjectQuantity
    if Sorted_DecisionValues(i) ~= DecisionValue_prev
        AUC = AUC + (FP - FP_prev) * ((TP + TP_prev) / 2);
        DecisionValue_prev = Sorted_DecisionValues(i);
        TP_prev = TP;
        FP_prev = FP;
        
        TP_Array = [TP_Array TP/P];
        FP_Array = [FP_Array FP/N];
        
        Accuracy_Array = [Accuracy_Array (TP + N - FP) / (P + N)];
    end
    if Label(OriginPos(i)) == -1
        TP = TP + 1;
    else
        FP = FP + 1;
    end 
end
AUC = AUC + (FP - FP_prev) * ((TP + TP_prev) / 2);
AUC = AUC / (length(find(Label == 1)) * length(find(Label == -1)));

LargestAccuracy = max(Accuracy_Array);

if ROC_Draw_Flag

    TP_Array = [TP_Array TP/P];
    FP_Array = [FP_Array FP/N];

    LargestAccuracyIndex = find(Accuracy_Array == max(Accuracy_Array));
    disp(max(Accuracy_Array));
    OptimalIndex = find(TP_Array == max(TP_Array(LargestAccuracyIndex)));
    OptimalIndex = intersect(OptimalIndex, LargestAccuracyIndex);
    % Create axes
    figure1 = figure;
    axes1 = axes('Parent', figure1,...
        'FontSize', 20, ...
        'LineWidth', 2,  'CameraViewAngle', 6.5);
    plot(FP_Array, TP_Array, '-r','LineWidth', 3);
    hXLabel = xlabel('1 - Specificity');
    hYLabel = ylabel('Sensitivity');

%     for i = 1:length(OptimalIndex)
%         hold on;
%         plot(FP_Array(OptimalIndex(i)), TP_Array(OptimalIndex(i)));
%              text(FP_Array(OptimalIndex(i)) + 0.035, TP_Array(OptimalIndex(i)) - 0.05, ['Accuracy=' ...
%             num2str(Accuracy_Array(OptimalIndex(i)) * 100) '%'], 'FontSize', 16, ...
%             'FontName', 'Arial', 'Color', [0 0 0]);
%         
%     end
    set(gca, 'XTick', [0:0.2:1]);
    set(gca, 'YTick', [0:0.2:1]);
    set(gca, 'FontName', 'Arial', 'FontSize', 16);
    set(gca, 'Box', 'off', ...                                       % 边框
         'LineWidth',2,...
         'XGrid', 'off', 'YGrid', 'off', ...                     % 网格
         'TickDir', 'in', 'TickLength', [.0 .0], ...          % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...           % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])    % 坐标轴颜色
    set([hXLabel,hYLabel], 'FontName',  'Arial', 'FontSize', 18,'FontWeight' , 'bold')
    fileout = 'ROC_Result_Red';
    print(figure1,[result_dir '\' fileout,'.png'],'-r300','-dpng');
end


