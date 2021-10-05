function [replaced, replacement, features, percImprovement] = nicheCompete(newInds,fitness,map,d,p,features)
%nicheCompete - results of competition with map's existing elites
%
% Syntax:  [replaced, replacement, features, percImprovement] = nicheCompete(newInds,fitness,phenotypes,map,d,p,varargin)
%
% Inputs:
%   newInds - [NXM]     - New population to compete for niches
%   fitness - [NX1]     - Fitness values fo new population
%   map     - struct    - Population archive
%   d       - struct    - Domain definition
%   p       - struct    - QD configuration
%
% Outputs:
%   replaced    - [NX1] - Linear index of map cells to recieve replacements
%   replacement - [NX1] - Index of newInds to replace current elites in niche
%   features    - [NxNumFeatures] - return features
%
%
% Other m-files required: getBestPerCell.m
%
% See also: createMap, getBestPerCell, updateMap

% Author: Adam Gaier, Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de, alexander.hagg@h-brs.de
% Jun 2016; Last revision: 15-Aug-2019

%------------- BEGIN CODE --------------
[bestIndex, bestBin] = getBestPerCell(fitness, d, map.edges, features);
mapLinIndx = sub2ind(p.featureResolution,bestBin(:,1),bestBin(:,2));
% Compare to already existing samples
improvement = ~(fitness (bestIndex) < map.fitness(mapLinIndx)); % comparisons to NaN are always false
improvement(isnan(fitness(bestIndex))) = false;
replacement = bestIndex (improvement);
replaced    = mapLinIndx(improvement);

percImprovement = 0; %TODO


%------------- END OF CODE --------------