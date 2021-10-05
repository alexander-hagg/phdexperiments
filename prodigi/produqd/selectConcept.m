function order = selectConcept(conceptSizes, criterion, varargin)
%SELECTCONCEPT Summary of this function goes here
%   Detailed explanation goes here

%p.selectCriterion.type      = 'size';           % 'nearest', 'size', variance', 'sparsity', 'map'
%p.selectCriterion.valType   = 'rank';           % 'rank', 'percentile'
%p.selectCriterion.value     = [1];              % [1 2 3] ('type' top three), [0.95 1.0] ('type' top 5%)
                                                % [-1 -2 -3] ('type' bottom three), [0.0 0.05] ('type' bottom 5%)
%p.selectCriterion.oneTime   = false;            % If true, selection is only applied once, 
                                                % after which selection is based on proximity to prototype in first selection

%% 
% Sort according to criterion type
disp(['Selection by ' criterion.type]);

if strcmp(criterion.type, 'size')
    [~,order] = sort(conceptSizes, 'descend'); 
elseif strcmp(criterion.type, 'nearest')
    if nargin < 4
        disp(['ERROR: when selecting by ' criterion.type ' criterion: not enough input arguments']);
    else
        oldConceptID = varargin{1};oldPrototypes = varargin{2};newPrototypes = varargin{3};
        %disp(['prot old: ' mat2str(size(oldPrototypes)) ', prot new: ' mat2str(size(newPrototypes))]);
        distances = pdist2(newPrototypes,oldPrototypes(oldConceptID,:));
        [~,order] = sort(distances, 'ascend'); 
    end
    
elseif strcmp(criterion.type, 'variance')
    disp(['ERROR: ' criterion.type ' criterion not implemented']);
    
elseif strcmp(criterion.type, 'sparsity')
    disp(['ERROR: ' criterion.type ' criterion not implemented']);   
    
elseif strcmp(criterion.type, 'map')
    disp(['ERROR: ' criterion.type ' criterion not implemented']);

elseif strcmp(criterion.type, 'prototypedist')
    if nargin < 2
        disp(['ERROR: when selecting by ' criterion.type ' criterion: not enough input arguments']);
    else
        prototypes = varargin{1};
        protDists = pdist2(prototypes,prototypes);
        [~,maxID] = max(protDists(:));
        [I,J] = ind2sub(size(protDists),maxID);
        order = [I J];
    end
end

if strcmp(criterion.valType, 'rank')
    order = order(criterion.value);    
elseif strcmp(criterion.valType, 'percentile')
    disp(['ERROR: ' criterion.valType ' criterion value type not implemented']);
end

end





