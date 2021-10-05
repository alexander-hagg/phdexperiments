function [children,selectedLinIDs] = createChildren(map, p, d)
%createChildren - produce new children through mutation of map elite
%
% Syntax:  children = createChildren(map,p)
%
% Inputs:
%   map - Population struct
%    .fitness
%    .genes
%    .<additional info> (e.g., drag, lift, etc)
%   p   - SAIL hyperparameter struct
%    .nChildren - number of children created
%    .mutSigma  - sigma of gaussian mutation applied to children
%   d   - Domain description struct
%    .dof       - Degrees of freedom (genome length)
%
% Outputs:
%   children - [nChildren X genomeLength] - new solutions
%
%

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 01-Aug-2017

%------------- BEGIN CODE --------------

% Remove empty bins from parent pool
parentPool = reshape(map.genes,[numel(map.fitness), d.dof]);
valid = ~isnan(parentPool(:,1));
parentPool(~valid,:) = [];

% Choose parents and create mutation
if strcmp(p.selectProcedure,'random')
    selection = randi([1 size(parentPool,1)], [p.nChildren 1]);
elseif strcmp(p.selectProcedure,'curiousness')
    validMapLinids = 1:length(valid);
    validMapLinids(~valid) = [];
    % Setup tournament selection
    selection = randi([1 size(parentPool,1)], [p.nChildren 2]);
    selectedLinIDs = validMapLinids(selection);
    
    curiousnessParents = reshape(map.curiousness,[numel(map.fitness), 1]);
    curiousnessParents(~valid) = [];
    curiousnessParents = curiousnessParents(selection);
    
    % Get most evolvable parents
    [~,sorted] = sort(curiousnessParents,2,'descend');
    idx = sub2ind(size(curiousnessParents), (1:size(sorted,1))', sorted(:,1));
    
    selectedLinIDs = selectedLinIDs(idx);
    selection = selection(idx);    
end
parents = parentPool(selection, :);

mutation = randn(p.nChildren,d.dof) .* p.mutSigma;

% Apply mutation
children = parents + mutation;
children(children>d.ranges(2)) = d.ranges(2); children(children<d.ranges(1)) = d.ranges(1);



%------------- END OF CODE --------------