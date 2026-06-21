function [resampled_structure, resampled_function,indices] = random_uniform_sample(structure_feature, varfeature,n_bins)
% RANDOM_UNIFORM_SAMPLE
% 在 data_vector 中进行“随机且均匀”的抽样，每次抽 n_bins 个点
%
% 输入：
%   structure_feature : [N×4] 结构特征矩阵
%   varfeature        : [N×1] 功能排序变量（用于决定结构排序）
%   n_bins      : 想要的采样数量（例如 190）

%
% 输出：
%   resampled_structure : [n_bins×4] 每 bin 中选出的结构特征
%   resampled_function  : [n_bins×1] 对应功能特征
if nargin < 3
    n_bins = 190;
end

n_total = length(varfeature);
assert(size(structure_feature,1) == n_total, '输入维度不匹配');


% 分 bin（均匀分隔）
bin_edges = round(linspace(1, n_total+1, n_bins+1));  % 生成 bin 边界

resampled_structure = zeros(n_bins, 4);
resampled_function = zeros(n_bins, 1);
indices = zeros(n_bins, 1);


for i = 1:n_bins
    idx_start = bin_edges(i);
    idx_end = bin_edges(i+1) - 1;

    % 防止 bin 太小
    if idx_start <= idx_end
        idx = randsample(idx_start:idx_end, 1);  % 在子区间内随机取一个
    else
        idx = idx_start;  % 只有一个点
    end

    indices(i) = idx;
    % 可选取策略：idx_start（第一个）、中间、随机
    resampled_structure(i,:) = structure_feature(idx, :);
    resampled_function(i) = varfeature(idx);
end


