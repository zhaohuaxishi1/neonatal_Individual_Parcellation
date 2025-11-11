%% visualize_structural_age_effects.m
% Visualize GAM-derived age effects on neonatal cortical structure
% Outputs: .shape.gii surface maps for uncorrected F statistics
% Author: jlzhao

clear; clc;

%% ---------------------------- Paths ----------------------------
age_effect_dir = fullfile('output', 'gam_structurematuration');
fsaverage_label_dir = './surface_template/fsaverage4/label';
template_dir = './surface_template/standard_mesh_atlases/resample_fsaverage';

% Load cortical label annotations
lh_annot = fullfile(fsaverage_label_dir, 'lh.aparc.a2009s.annot');
rh_annot = fullfile(fsaverage_label_dir, 'rh.aparc.a2009s.annot');
[~, lh_labels, ~] = read_annotation(lh_annot);
[~, rh_labels, ~] = read_annotation(rh_annot);
all_labels = [lh_labels; rh_labels];  % Combined left and right hemisphere labels

hemi_names = {'left', 'right'};
hemi_short = {'L', 'R'};

%% ---------------------------- Uncorrected Effects ----------------------------
fprintf('===> Processing uncorrected absolute F-statistics...\n');

metric_names_uncorrected = {'SulcalDepth','Curvature','Myelination','Thickness'};
n_metrics_uncorrected = length(metric_names_uncorrected);
n_vertices = 5124;

age_f_matrix = zeros(n_metrics_uncorrected, n_vertices);

% Create output folder
uncorr_out_dir = fullfile(age_effect_dir, 'UnthreshAbs');
if ~exist(uncorr_out_dir, 'dir'), mkdir(uncorr_out_dir); end

for i = 1:n_metrics_uncorrected
    data = load(fullfile(age_effect_dir, ['AgeEffect_Structure_', metric_names_uncorrected{i}, '.mat']));
    f_vals = data.Gam_F_Vector_All;
    f_vals(all_labels == 1644825) = 0;  % Mask medial wall
    age_f_matrix(i, :) = f_vals;
end

age_f_abs = abs(age_f_matrix)';  % [vertices x metrics]
save(fullfile(uncorr_out_dir, 'Gam_Age_Abs.mat'), 'age_f_abs');

% Surface visualization (unthresholded)
for m = 1:n_metrics_uncorrected
    for h = 1:2
        template_file = fullfile(template_dir, ...
            sprintf('fsaverage4.%s.midthickness_va_avg.3k_fsavg_%s.shape.gii', hemi_short{h}, hemi_short{h}));
        output_file = fullfile(uncorr_out_dir, ...
            sprintf('AgeEffect_%s_hemi-%s_mesh-fsaverage4.shape.gii', metric_names_uncorrected{m}, hemi_names{h}));

        g = gifti(template_file);
        if h == 1
            g.cdata = age_f_abs(1:2562, m);
        else
            g.cdata = age_f_abs(2563:end, m);
        end
        save(g, output_file);
    end
end

