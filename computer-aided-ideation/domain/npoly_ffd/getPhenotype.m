function phenotypes = getPhenotype(genomes,d)
%getPhenotype - Express one or more genomes
%
% Syntax:  phenotypes = getPhenotype(genomes,d)
%
% Inputs:
%    genomes        - [NxM] - N genomes with dof = M
%    d              - struct - Domain description struct
%
% Outputs:
%    phenotypes          - cell[Nx1] - Cell array containing all phenotypes
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Aug 2019; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

for x=1:size(genomes,1)
    genome = genomes(x,:);
    if ~isnan(genome(1))
        theta = d.base(1,:) + 0.5*genome(1:size(d.base,2));
        rho = d.base(2,:) + genome(size(d.base,2)+1:end);
        [xCoords,yCoords] = pol2cart(theta,rho);
        phenotypes{x} = polyshape(xCoords,yCoords);
    else
        phenotypes{x} = nan;
    end
end
end

%------------- END CODE --------------