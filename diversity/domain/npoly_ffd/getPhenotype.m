function phenotypes = getPhenotype(genomes,base)
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
tol = 0.01;

nGenomes = 1;
if size(genomes,2)>1
    nGenomes = size(genomes,1);
end

for x=1:nGenomes
    if size(genomes,2)>1
        genome = genomes(x,:);
    else
        genome = genomes';
    end
    if ~isnan(genome(1))
        x1 = base(1:end/2) + genome(1:end/2);
        y1 = base(end/2+1:end) + genome(end/2+1:end);
        phenotypes{x} = polyshape(x1,y1);
        %phenotypes{x} = rmslivers(phenotypes{x},tol);
    else
        phenotypes{x} = nan;
    end
end
end

%------------- END CODE --------------