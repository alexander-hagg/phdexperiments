function [replaced, replacement, percImprovement] = nicheCompete(newInds,fitness,values,map,d)
%nicheCompete - results of competition with map's existing elites
%
% Syntax:  [replaced, replacement] = nicheCompete(newInds,fitness,map,p)
%
% Inputs:
%   newInds - [NXM]     - New population to compete for niches
%   fitness - [NX1]     - Fitness values fo new population
%   map     - struct    - Population archive
%   d       - struct    - Domain definition
%
% Outputs:
%   replaced    - [NX1] - Linear index of map cells to recieve replacements
%   replacement - [NX1] - Index of newInds to replace current elites in niche
%
% Example:
%
% Other m-files required: getBestPerCell.m
%
% See also: createMap, getBestPerCell, updateMap

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 02-Aug-2017

%------------- BEGIN CODE --------------
[bestIndex, bestBin] = getBestPerCell(newInds,fitness, values, d, map.edges);
mapLinIndx = sub2ind(d.featureRes,bestBin(:,1),bestBin(:,2));

% Compare to already existing samples
improvement = ~(fitness (bestIndex) >= map.fitness(mapLinIndx)); % comparisons to NaN are always false
improvement(isnan(fitness(bestIndex))) = false;
replacement = bestIndex (improvement);
replaced    = mapLinIndx(improvement);
%------------- END OF CODE --------------