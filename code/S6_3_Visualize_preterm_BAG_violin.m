%% ------------------------------------------------------------------------
% visualize_preterm_BAG_violin.m
%
% Visualization of personalized brain age gap (BAG) differences between
% term-born and preterm infants.
%
% This script generates violin plots for:
%   1. Whole-brain BAG
%   2. Network-level BAG
%
% Networks shown:
%   - Whole brain
%   - DMN-p
%   - DAN
%   - Superior Parietal
%   - Frontal-Limbic
%
% Corresponding manuscript:
% "Personalized Functional Network Topographic Mapping in Human Neonates"
%
% Author: jlzhao
%% ------------------------------------------------------------------------

clc;
clear;

%% --------------------------- Configuration ---------------------------

data_dir   = fullfile('data','brain_age_BAG');
output_dir = fullfile('output','Figure7_BAG');

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% --------------------------- Load Data ---------------------------

load(fullfile(data_dir,'Network_Level_Gap.mat'));

preterm_bag_network = gap_resid_preterm;
term_bag_network    = gap_resid_term;

% Network ordering used throughout the manuscript
network_order = [11 10 4 5 2 9 8 7 6 3 1];

preterm_bag_network = ...
    preterm_bag_network(:,network_order);

term_bag_network = ...
    term_bag_network(:,network_order);

% Selected significant networks
selected_networks = [3 4 6 11];

preterm_bag_network = ...
    preterm_bag_network(:,selected_networks);

term_bag_network = ...
    term_bag_network(:,selected_networks);

%% --------------------------- Whole-Brain BAG ---------------------------

load(fullfile(data_dir,'gap_resid_preterm_term.mat'));

preterm_bag = ...
    [gap_resid_preterm, preterm_bag_network];

term_bag = ...
    [gap_resid_term, term_bag_network];

clear gap_resid_preterm gap_resid_term

%% --------------------------- Color Settings ---------------------------

color_map = [
    0 100 120
    120 180 190
    10  83  73
    16 139 122
    79 113  75
    132 188 125
    136  36  30
    227  60  51
    136 100  45
    228 167  76
    ] ./ 255;

%% --------------------------- Figure Settings ---------------------------

figure_width  = 10;
figure_height = 6;

figure('Units','centimeters',...
    'Position',[0 0 figure_width figure_height]);

hold on

%% --------------------------- Violin Plot ---------------------------

for network_id = 1:5
    
    preterm_values = preterm_bag(:,network_id);
    term_values    = term_bag(:,network_id);
    
    x_preterm = 0.5 + 3*(network_id-1);
    x_term    = 1.4 + 3*(network_id-1);
    
    Violin({preterm_values}, x_preterm,...
        'ViolinColor',{color_map(2*(network_id-1)+2,:)},...
        'HalfViolin','full',...
        'ShowData',false,...
        'ShowBox',true,...
        'ShowMedian',true,...
        'ViolinAlpha',{0.9});
    
    Violin({term_values}, x_term,...
        'ViolinColor',{color_map(2*(network_id-1)+2,:)},...
        'HalfViolin','full',...
        'ShowData',false,...
        'ShowBox',true,...
        'ShowMedian',true,...
        'ViolinAlpha',{0.9});
    
    
end

%% --------------------------- Axis Settings ---------------------------

for network_id = 1:5
    
    x_preterm(network_id) = ...
        0.5 + 3*(network_id-1);
    
    x_term(network_id) = ...
        1.4 + 3*(network_id-1);
    
end

set(gca,...
    'Box','off',...
    'LineWidth',1,...
    'TickDir','out',...
    'TickLength',[.005 .005],...
    'FontName','Arial',...
    'FontSize',14);

set(gca,...
    'XTick',(x_preterm+x_term)/2,...
    'XTickLabel',{'All','3','4','6','11'},...
    'XLim',[-0.5 15],...
    'YTick',-20:10:20,...
    'YLim',[-20 20]);

ylabel('Brain Age Gap (weeks)',...
    'FontSize',14,...
    'FontName','Arial');

%% --------------------------- Save Figure ---------------------------

set(gcf,...
    'PaperUnits','centimeters',...
    'PaperPosition',[0 0 figure_width figure_height]);

print(gcf,...
    fullfile(output_dir,...
    'Figure7_BAG_Violin.png'),...
    '-dpng','-r300');

fprintf('Figure saved successfully.\n');
