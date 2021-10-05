function [fitness,polygons,rawFitness,relativeSelectionDistance] = fitfun(genomes,d)
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


rawFitness = zeros(length(polygons),1);

for i=1:length(polygons)
    %disp([int2str(i) '/' int2str(length(polygons))]);
    % Return NaN for invalid phenotypes
    if isa(polygons{i},'double') || polygons{i}.NumRegions == 0
        rawFitness(i,1) = nan;
    else
        % Point Symmetry
        pt = polygons{i}.Vertices;
        pt = pt - mean(pt); % Center
        a = pt(1:ceil(end/2),:); % Take first half of points
        x = pt(ceil(end/2)+mod(size(pt,1)+1,2):end,:); % Second half of points
        rawFitness(i) = mean(sqrt(sum(((a+x)'.^2),1)));
        
    end
end

symmetryFitness = 1./(1 + rawFitness);


% Involve the user selection only if it exists as a field in the d struct
if ~isfield(d,'selection') || ~isfield(d.selection,'models') || isempty(d.selection.models)
    fitness = symmetryFitness;
    rawFitness = symmetryFitness;
else
    fitness = symmetryFitness;
    for i=1:length(d.selection.models)
        % Get features
        hypersurfaceID = i;
        features = predictFeatures(polygons,d.selection.models{hypersurfaceID});
        % Get phenotypic distances to user selection
        phenoDistToSelection = pdist2(features,d.selection.selected{hypersurfaceID});
        phenoDistToDeselection = pdist2(features,d.selection.deselected{hypersurfaceID});
        
        % Get total distance to entire selection
        % Prefer solutions that are "in between" solutions rather than the ones
        % only close to one of the selected solutions
        selectionDistances = min(phenoDistToSelection,[],2);
        deselectionDistances = min(phenoDistToDeselection,[],2);
        
        relativeSelectionDistance = selectionDistances./(selectionDistances+deselectionDistances);
        %penalty = relativeSelectionDistance;
        % Violation: relativeSelectionDistance > 0.5
        %penalty(penalty < 0.5) = 0.5;
        %selectionFitness = 1-(penalty - 0.5)*2;
        %penalty(penalty>0.5) = 1;
        selectionFitness = 1 - relativeSelectionDistance;
        
        fitness = fitness .* selectionFitness;
    end
    rawFitness = [fitness,selectionFitness];
end

% Limit fitness between 0 and 1
fitness(fitness>1) = 1;
fitness(fitness<0) = 0;

end

%------------- END CODE --------------