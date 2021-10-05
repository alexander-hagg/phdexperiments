function [replaced, replacement, features, percImprovement] = nicheCompete(newInds,fitness,map,d,p,features)
%nicheCompete - results of competition with map's existing elites
%
% Syntax:  [replaced, replacement, features, percImprovement] = nicheCompete(newInds,fitness,map,d,p,varargin)
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

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jul 2019; Last revision: 15-Aug-2019

%------------- BEGIN CODE --------------
replacement = false(size(newInds,1),1);

percImprovement = 0;


% Get distance to elites
if ~isempty(map.genes)
    eliteDistance = pdist2(features,map.features);
else
    eliteDistance = [];
end

% Get distance between candidates
distances = [eliteDistance pdist2(features,features)];
distances(distances==0) = nan; %TODO: this is a hack to prevent comparisons of a candidate with itself

% Compete if needed
competing = distances < map.config.competeDistance;%.*fitnessScaling;
competition = ([map.fitness, fitness] .* competing);
competition(~competing) = nan;

% Add competing candidates that improve the map
won = fitness' > competition;
takehome = won;
takehome(~competing) = 1; % Add non-competing
replacement = all(takehome'==1);

replaced = false(length(map.fitness),1);
if ~isempty(map.genes)
    % only remove elites from the map here
    distanceCompetition = (distances .* competing);
    distanceCompetition(~competing) = nan;
    
    % Get all distances
    % Only distances that were won by new candidates
    removalCandidates = distanceCompetition.*won;
    % Ignore non-competitions
    removalCandidates(removalCandidates==0) = nan;
    % Get only candidate distances from map
    removalCandidates = removalCandidates(:,1:end-length(fitness));
    [~, nn] = min(removalCandidates');
    removeIDs = nn(~all(isnan(removalCandidates')));
    if ~isempty(removeIDs); replaced(removeIDs) = 1;end
end

%% TODO MAX BINS
%p.maxBins
%------------- END OF CODE --------------