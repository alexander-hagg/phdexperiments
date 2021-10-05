function [fitness,polygons,rawFitness] = fitfun(genomes,d)
%fitfun - pointSymmetry fitness function
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

if isempty(genomes)
    fitness = [];
    polygons = [];
    rawFitness = [];
    return;
end


[polygons] = d.getPhenotype(genomes);
[~,boolmaps] = getPhenotypeBoolean(polygons,d.resolution);
rawFitness = zeros(length(polygons),1);

%%
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
    fitness(i) = rawFitness(i)./(maxFitness(i) + 1e-5);
%     figure(1);hold off;
%     imagesc(boolmaps{i})
%     hold on;
%     scatter(c(1),c(2),'g','filled')
%     scatter(BB.BoundingBox(1),BB.BoundingBox(2),32,'k','filled');
%     scatter(BB.BoundingBox(1)+BB.BoundingBox(3),BB.BoundingBox(2),32,'k','filled');
%     scatter(BB.BoundingBox(1),BB.BoundingBox(2)+BB.BoundingBox(4),32,'k','filled');
%     scatter(BB.BoundingBox(1)+BB.BoundingBox(3),BB.BoundingBox(2)+BB.BoundingBox(4),32,'k','filled');
%     scatter(cBB(1),cBB(2),32,'r','filled')
%     axis equal;
%     title(pdist2(c,cBB))
%     pause(1)
end

fitness = fitness';
% Limit fitness between 0 and 1
fitness(fitness>1) = 1;
fitness(fitness<0) = 0;

end

%------------- END CODE --------------