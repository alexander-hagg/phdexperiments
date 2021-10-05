function [fitness,phenotypes] = fitnessMaxAreaMinCircumference(genomes,d)
%npolyObjective - domain's single-objective fitness function
% Point symmetry
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
numSamplePts = 100;

[polygons,phenotypes] = d.getPhenotype(genomes);
fitness = zeros(length(polygons),1);

for i=1:length(polygons)
    % Return NaN for invalid phenotypes
    if isa(polygons{i},'double') || polygons{i}.NumRegions == 0
        fitness(i) = nan;
    else
        vertices = polygons{i}.Vertices;
        vertices = vertices(~all(isnan(vertices)'),:);
        vertices(end+1,:) = vertices(1,:);
        pt = interparc(numSamplePts+1,vertices(:,1),vertices(:,2),'linear');
        pt(end,:) = [];
        
        pol = polyshape(pt);
        areaDivCirc = area(pol)/perimeter(pol);
        [x0,y0] = centroid(pol);
        meanRadius = mean(pdist2(pol.Vertices(~all(isnan(pol.Vertices)'),:),[x0,y0]));
        fitness(i) = 2*areaDivCirc/meanRadius;

    end
end

fitness(fitness<0) = 0;
fitness(fitness>1) = 1;

fitness = fitness';

end

%------------- END CODE --------------