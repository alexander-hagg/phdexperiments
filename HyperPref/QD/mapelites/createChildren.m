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
parentPoolFitness = reshape(map.fitness,[],1);
valid = ~isnan(parentPool(:,1));
parentPool(~valid,:) = [];
parentPoolFitness(~valid) = [];

% Choose parents 
selection = randi([1 size(parentPool,1)], [p.nChildren 1]);
parents = parentPool(selection, :);
parentFitness = parentPoolFitness(selection);

if p.crossover
    children = [];
    crossoverChance = 0.6;
    for i=1:size(parents,1)/2
        pair = parents([i,i+size(parents,1)/2],:);
        if parentFitness(i)< parentFitness(i+size(parents,1)/2)
            pair = flipud(pair);
        end
        for j=1:2 % 2 children per parent pair
            crossoverRoll = rand(1,size(pair,2)) < crossoverChance;
            children(end+1,:) = crossoverRoll.*pair(1,:) + ~crossoverRoll.*pair(2,:);
        end
    end
else
    children = parents;
end

% Mutation
mutation = randn(p.nChildren,d.dof) .* p.mutSigma .* range(d.ranges');
children = children + mutation;

[id] = find(children>d.ranges(:,2)');
[~,col] = find(children>d.ranges(:,2)');
children(id) = d.ranges(col,2);

[id] = find(children<d.ranges(:,1)');
[~,col] = find(children<d.ranges(:,1)');
children(id) = d.ranges(col,2);


%------------- END OF CODE --------------