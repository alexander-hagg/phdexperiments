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
[~,boolmaps] = getPhenotypeBoolean(polygons,d.resolution);

rawFitness = zeros(length(polygons),1);

for i=1:length(boolmaps)
    c = regionprops(boolmaps{i},'centroid');
    c = c.Centroid;
    BB = regionprops(boolmaps{i},'BoundingBox');
    maxID = 1; maxSize = BB(1).BoundingBox(3)*BB(1).BoundingBox(4);
    if numel(BB) > 1
        for j=2:numel(BB)
            tsize = BB(j).BoundingBox(3)*BB(j).BoundingBox(4);
            if tsize > maxSize
                maxID = j; maxSize = tsize;
            end
        end
    end
    BB = BB(maxID);
    cBB = [BB.BoundingBox([1,2]) + 0.5 * BB.BoundingBox([3,4])];
    rawFitness(i) = pdist2(c,cBB);
    maxFitness(i) = sqrt(BB.BoundingBox(3).^2+BB.BoundingBox(4).^2);
    asymmetryFitness(i) = rawFitness(i)./(maxFitness(i) + 1e-5);
end

% Get phenotypic distances to user selection
phenoDistToSelection = pdist2(features,d.selectedShapes);
phenoDistToDeselection = pdist2(features,d.deselectedShapes);

% Get total distance to entire selection
% Prefer solutions that are "in between" solutions rather than the ones
% only close to one of the selected solutions
selectionDistances = min(phenoDistToSelection,[],2);
deselectionDistances = min(phenoDistToDeselection,[],2);


relativeSelectionDistance = selectionDistances./(selectionDistances+deselectionDistances);
if strcmp(d.selectPenalty,'none')
    % No constraints
    selectionFitness = ones(size(relativeSelectionDistance,1),1);
elseif strcmp(d.selectPenalty,'inverseDistance')
    % Calculate fitness from selection distances
    selectionFitness = 1./(1+selectionDistances);
elseif strcmp(d.selectPenalty,'relativeDistance')
    selectionFitness = 1-relativeSelectionDistance;
elseif strcmp(d.selectPenalty,'relativeDistanceOnlyPenalizeConstraintViolation')
    penalty = relativeSelectionDistance;
    % Violation: relativeSelectionDistance > 0.5
    penalty(penalty < 0.5) = 0.5;    
    selectionFitness = 1-(penalty - 0.5)*2;
end

asymmetryFitness = asymmetryFitness';

fitness = asymmetryFitness .* selectionFitness;

% Limit fitness between 0 and 1
fitness(fitness>1) = 1;
fitness(fitness<0) = 0;

rawFitness = [asymmetryFitness,selectionFitness];

end

%------------- END CODE --------------