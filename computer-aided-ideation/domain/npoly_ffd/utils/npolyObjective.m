function [fitness,values,phenotypes] = npolyObjective(genomes,d)
%npolyObjective - domain's single-objective fitness function
% The GUI currently works with normalized fitness only.
%
% Syntax:  [fitness,values,phenotypes] = npolyObjective(genomes,d)
%
% Inputs:
%    genomes        - [NxM] - N genomes with dof = M
%    d              - cell - Domain configuration.
%
% Outputs:
%    fitness        - [Nx1] - Validation flags
%    values         - cell[?] - may contain extra values extracted from
%                               fitness function
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

values = [];
phenotypes = getPhenotype(genomes,d);
fitness = zeros(length(phenotypes),1);
for i=1:length(phenotypes)
    if phenotypes{i}.NumRegions == 0
        fitness(i) = nan;
    else
        % Radial Symmetry
        [xc,yc] = phenotypes{i}.centroid;
        vertices = phenotypes{i}.Vertices;
        vertices = vertices(~all(isnan(vertices)'),:);
        vertices = vertices -[xc,yc];
        a = vertices(1:end/2,:);
        x = -1*(vertices(end/2+1:end,:));
        distances = hypot(a(:,1)-x(:,1),a(:,2)-x(:,2));
        fitness(i) = sum(distances);
    end
end

dof = size(genomes,2);
fitness = 1-(fitness'./(5*dof/8));  % try to normalize between 0 and 1, values depend on DOF
                                    % (approximately correct)

fitness(fitness<0) = 0;
fitness(fitness>1) = 1;

end

%------------- END CODE --------------