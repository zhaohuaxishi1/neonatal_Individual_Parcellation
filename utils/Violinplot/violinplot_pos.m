function violins = violinplot(data, cats, varargin) 
%Violinplots plots violin plots of some data and categories
%   [文档略，保留原注释] 

% Copyright (c) 2016, Bastian Bechtold
% This code is released under the terms of the BSD 3-clause license

hascategories = exist('cats','var') && not(isempty(cats));

% Parse the optional grouporder argument
grouporder = {};
idx=find(strcmp(varargin, 'GroupOrder'));
if ~isempty(idx) && numel(varargin)>idx
    if iscell(varargin{idx+1})
        grouporder = varargin{idx+1};
        varargin(idx:idx+1)=[];
    else
        error('Second argument of ''GroupOrder'' optional arg must be a cell of category names')
    end
end

% Parse optional X argument (newly added)
xpositions = [];
idx = find(strcmp(varargin, 'X'));
if ~isempty(idx) && numel(varargin)>idx
    xpositions = varargin{idx+1};
    varargin(idx:idx+1) = [];
end

% ViolinColor/Alpha 参数检查略...

% tabular data
if isa(data, 'dataset') || isstruct(data) || istable(data)
    % [略]
    for n=1:length(catnames)
        thisData = data.(catnames{n});
        xpos = n;
        if ~isempty(xpositions)
            xpos = xpositions(n);
        end
        violins(n) = Violin({thisData}, xpos, varargin{:});
    end
    set(gca, 'XTick', 1:length(catnames), 'XTickLabels', catnames);
    set(gca,'Box','on');
    return
elseif iscell(data) && length(data(:))==2
    % [略]
elseif iscell(data) && length(data(:))>2
    error('Up to two datasets can be compared');
elseif isnumeric(data)
    if hascategories && numel(data) == numel(cats)
        if isempty(grouporder)
            cats = categorical(cats);
        else
            cats = categorical(cats, grouporder);
        end

        catnames = (unique(cats));
        catnames_labels = {};
        for n = 1:length(catnames)
            thisCat = catnames(n);
            catnames_labels{n} = char(thisCat);
            thisData = data(cats == thisCat);
            xpos = n;
            if ~isempty(xpositions)
                xpos = xpositions(n);
            end
            violins(n) = Violin({thisData}, xpos, varargin{:});
        end
        set(gca, 'XTick', 1:length(catnames), 'XTickLabels', catnames_labels);
        set(gca,'Box','on');
        return
    else
        data = {data};
    end
end

% 1D data, no categories
if ~hascategories && isvector(data{1})
    xpos = 1;
    if ~isempty(xpositions)
        xpos = xpositions(1);
    end
    violins = Violin(data, xpos, varargin{:});
    set(gca, 'XTick', xpos);
% 2D data with or without categories
elseif ismatrix(data{1})
    for n=1:size(data{1}, 2)
        thisData = cellfun(@(x)x(:,n),data,'UniformOutput',false);
        xpos = n;
        if ~isempty(xpositions)
            xpos = xpositions(n);
        end
        violins(n) = Violin(thisData, xpos, varargin{:});
    end
    set(gca, 'XTick', 1:size(data{1}, 2));
    if hascategories && length(cats) == size(data{1}, 2)
        set(gca, 'XTickLabels', cats);
    end
end

set(gca,'Box','on');

end
