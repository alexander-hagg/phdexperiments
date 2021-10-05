function [fitness,phenotypes] = fitnessPointSymmetry(genomes,d)
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
numSamplePts = 10;
        
[phenotypes,polygons] = d.getPhenotype(genomes);
fitness = zeros(length(polygons),1);

for i=1:length(polygons)
    % Return NaN for invalid phenotypes
    if isa(polygons{i},'double') || polygons{i}.NumRegions == 0
        fitness(i,1) = nan;
    else
        % Point Symmetry
        vertices = polygons{i}.Vertices;
        vertices = vertices(~all(isnan(vertices)'),:);
        vertices(end+1,:) = vertices(1,:);
        pt = interparc(numSamplePts+1,vertices(:,1),vertices(:,2),'linear');
        pt(end,:) = [];
        allVertices(i,:,:) = pt;
        a = pt(1:end/2,:);
        x = pt(end/2+1:end,:);
        
        fitness(i) = sum(sqrt(sum((a+x).^2,2)));
    end
end

fitness = 1./(1+fitness');
fitness(fitness<0) = 0;
fitness(fitness>1) = 1;

end

%------------- END CODE --------------