function children = createChildren(map, p, d)
%createChildren - produce new children through mutation of elites from map.
%
% Syntax:  children = createChildren(map,p)
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

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jul 2019; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

% Remove empty bins from parent pool
parentPool = map.genes;

% Choose parents and create mutation
selection = randi([1 size(parentPool,1)], [p.nChildren 1]);
parents = parentPool(selection, :);

% Apply mutation
mutation = randn(p.nChildren,d.dof) .* p.mutSigma;
children = parents + mutation;
children(children>d.ranges(2)) = d.ranges(2); children(children<d.ranges(1)) = d.ranges(1);



%------------- END OF CODE --------------