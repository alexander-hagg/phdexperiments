function [fitness,features,polygons,rawFitness] = fitfun(genomes,d)
%fitfun - "ui circle or triangle" fitness function
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
%% ------------- BEGIN CODE --------------

if isempty(genomes)
    fitness = [];
    polygons = [];
    rawFitness = [];
    return;
end


[polygons] = d.getPhenotype(genomes);
rawFitness = zeros(length(polygons),1);

for i=1:length(polygons)
    % Return NaN for invalid phenotypes
    if isa(polygons{i},'double') || polygons{i}.NumRegions == 0
        rawFitness(i,1) = nan;
    else
        pt = polygons{i}.Vertices;
        pt = pt - mean(pt); % Center
        distances = sqrt(pt(:,1).^2+pt(:,2).^2); 
        
        % Circleness
        circleness(i) = std(distances); % Variance of radii
        
        % Triangleness
        % mean distance (to be maximized)
        [peaks,locs,w,p] = findpeaks([distances(end-5:end);distances]);
        isTriangle = ~logical(abs(numel(peaks)-3)); % Is a triangle
        
        triangleness(i) = isTriangle*mean(p)./max(distances(:)); % Sum of Peak Prominences
    end
end

% Add two fitness functions, but do not rescale between 0 and 1. Instead,
% cap the values between 0 and 1 to allow full triangles and circles
% (ignoring the other fitness function)

circleness = 2./(1+5*circleness)-1;
circleness(circleness>1) = 1;
circleness(circleness<0) = 0;
triangleness(triangleness>1) = 1;
triangleness(triangleness<0) = 0;

fitness = (circleness) + triangleness;
rawFitness = [circleness',triangleness'];
fitness = fitness';
fitness(fitness>1) = 1;
fitness(fitness<0) = 0;

features = categorize(polygons,d);

end

%------------- END CODE --------------