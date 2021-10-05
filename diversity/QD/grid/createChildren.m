function [children,selectedLinIDs] = createChildren(map, p, d)
%createChildren - produce new children through mutation of elites from map.
% Selection: randomly taken from the map
%
% Syntax:  children = createChildren(map, p, d)
%
% Inputs:
%   map - struct - Population struct
%      .fitness
%      .genes
%      .<additional info> (e.g., drag, lift, etc)
%   p   - struct - QD configuration struct
%      .nChildren - number of children created
%      .mutSigma  - sigma of gaussian mutation applied to children
%   d   - struct - Domain description struct
%    .dof       - Degrees of freedom (genome length)
%
% Outputs:
%   children - [p.nChildren X d.dof] - new solutions
%
%

% Author: Adam Gaier, Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de, alexander.hagg@h-brs.de
% Jun 2016; Last revision: 15-Aug-2019

%------------- BEGIN CODE --------------

% Remove empty bins from parent pool
parentPool = reshape(map.genes,[numel(map.fitness), d.dof]);
valid = ~isnan(parentPool(:,1));
parentPool(~valid,:) = [];

% Choose parents and create mutation
if strcmp(p.selectProcedure,'random')
    selection = randi([1 size(parentPool,1)], [p.nChildren 1]);
elseif strcmp(p.selectProcedure,'curiosity')
    validMapLinids = 1:length(valid);
    validMapLinids(~valid) = [];
    % Setup tournament selection
    selection = randi([1 size(parentPool,1)], [p.nChildren 2]);
    selectedLinIDs = validMapLinids(selection);
    
    curiousParents = reshape(map.curiosity,[numel(map.fitness), 1]);
    curiousParents(~valid) = [];
    curiousParents = curiousParents(selection);
    
    % Get most evolvable parents
    [~,sorted] = sort(curiousParents,2,'descend');
    idx = sub2ind(size(curiousParents), (1:size(sorted,1))', sorted(:,1));
    
    selectedLinIDs = selectedLinIDs(idx);
    selection = selection(idx);    
end
parents = parentPool(selection, :);



% Mutation (no crossover!)
mutation = randn(p.nChildren,d.dof) .* p.mutSigma .* range(d.ranges');
children = parents + mutation;

[id] = find(children>d.ranges(:,2)');
[~,col] = find(children>d.ranges(:,2)');
children(id) = d.ranges(col,2);

[id] = find(children<d.ranges(:,1)');
[~,col] = find(children<d.ranges(:,1)');
children(id) = d.ranges(col,2);


%------------- END OF CODE --------------