%% ------------------------------------------------------------------------
% Network-level PLS Correlation Visualization
%
% Author: Jianlong Zhao
%
% This script visualizes the distribution of bootstrap-derived
% PLS correlations across functional networks using violin plots.
%
% Input:
%   str_func_corr_value.mat
%
% Output:
%   Step_2_network_level_PLS_All_Age.png
%% ------------------------------------------------------------------------

clc;
clear;

%% --------------------------- Load Data ---------------------------
% Paths
data_dir = 'data';

load(fullfile(data_dir, 'str_func_corr_value.mat'));   % r_value: bootstrap PLS correlations

result_dir = fullfile('output', 'PLS_Struc_Func_Age');
if ~exist(result_dir, 'dir'); mkdir(result_dir); end

%% --------------------------- Color Definition ---------------------------

C1 = [162, 20, 47];
C2 = [237, 177, 32];

network_colors = [
    120 180 190;
    C1;
    C2;
    C2;
    C2;
    C2;
    C1;
    C1;
    C2;
    C2;
    C2;
    C1];

network_colors = network_colors ./ 255;

%% --------------------------- Figure Settings ---------------------------

figureWidth  = 15;  % cm
figureHeight = 6;   % cm

figureHandle = figure;
set(figureHandle,...
    'Units','centimeters',...
    'Position',[0 0 figureWidth figureHeight]);

hold on;

%% --------------------------- Violin Plot ---------------------------

nNetworks = 12;

x_pos = 0.5 + 2*(0:nNetworks-1);

for net = 1:nNetworks

    network_values = r_value(:,net);

    Violin({network_values},...
        x_pos(net),...
        'ViolinColor',{network_colors(net,:)},...
        'HalfViolin','full',...
        'ShowData',false,...
        'ShowBox',true,...
        'ShowMedian',true,...
        'MarkerSize',2,...
        'ViolinAlpha',{0.9});

end

%% --------------------------- Axis Formatting ---------------------------

set(gca,...
    'Box','off',...
    'LineWidth',1,...
    'XGrid','off',...
    'YGrid','off',...
    'TickDir','out',...
    'TickLength',[0.005 0.005],...
    'XMinorTick','off',...
    'YMinorTick','off',...
    'XColor',[0.1 0.1 0.1],...
    'YColor',[0.1 0.1 0.1],...
    'FontName','Arial',...
    'FontSize',14);

set(gca,...
    'XTick',x_pos,...
    'XTickLabel',{'All','3','6','5','11','7','2','1','8','10','9','4'},...
    'XLim',[-0.5 24],...
    'YTick',0:0.2:0.8,...
    'YLim',[-0.05 0.8]);

ylabel('PLS correlation','FontSize',16,'FontName','Arial');

%% --------------------------- Export Figure ---------------------------

set(figureHandle,...
    'PaperUnits','centimeters',...
    'PaperPosition',[0 0 figureWidth figureHeight]);

output_file = [result_dir '\network_level_PLS_All_Age'];

print(figureHandle,...
    [output_file,'.png'],...
    '-dpng',...
    '-r300');

fprintf('Figure saved: %s.png\n',output_file);