function [fitness,polygons,rawFitness] = fitfun(genomes,d)
%fitfun - "ui compare to user selected shapes" fitness function
% Fitness is normalized between 0 and 1
%
% Syntax:  [fitness,phenotypes] = npolyObjective(genomes,d)
%
% Inputs:
%    genomes        - [NxM] - N genomes with dof = M
%    d              - cell - Domain configuration.
%
% Outputs:
%    fitness        - [Nx1] - Validation flags
%    phenotypes     - cell[Nx1] - phenotypes (to prevent recalculating
%                                 of phenotypes, we offer them back here
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jul 2019; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

if isempty(genomes)
    fitness = [];
    polygons = [];
    rawFitness = [];
    return;
end


[polygons] = d.getPhenotype(genomes);
% Get features
features = predictFeatures(polygons,d.userModel);

rawFitness = zeros(length(polygons),1);

% Get phenotypic distances to user selection
phenotypicDistances = pdist2(features,d.selectedShapes);

% Get total distance to entire selection
% Prefer solutions that are "in between" solutions rather than the ones
% only close to one of the selected solutions
selectionDistances = min(phenotypicDistances,[],2);

% Calculate fitness from selection distances
fitness = 1./(1+selectionDistances);

% Limit fitness between 0 and 1
fitness(fitness>1) = 1;
fitness(fitness<0) = 0;



end

%------------- END CODE --------------