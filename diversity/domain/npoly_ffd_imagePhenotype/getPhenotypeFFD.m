function polyshapes = getPhenotypeFFD(genomes,base)
%function [polyshapes,booleanMap,pixelCoordinates,flatbooleanMap] = getPhenotypeFFD(genomes,base,varargin)
%getPhenotype - Express one or more genomes
%
% Syntax:  polyshapes = getPhenotype(genomes,d)
%
% Inputs:
%    genomes        - [NxM] - N genomes with dof = M
%    d              - struct - Domain description struct
%
% Outputs:
%    polyshapes          - cell[Nx1] - Cell array containing all polyshapes
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Aug 2019; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------
nGenomes = 1; if size(genomes,2)>1; nGenomes = size(genomes,1); end
numSamplePts = 60;

for i=1:nGenomes
    if size(genomes,2)>1; genome = genomes(i,:); else;genome = genomes';end
    if ~isnan(genome(1))
        [theta,rho] = cart2pol(base(1:end/2)',base(end/2+1:end)');
        rho         = rho .* (genome(1:end/2)');
        theta       = theta + genome(end/2+1:end)';
        [x,y]       = pol2cart(repmat(theta,1,size(rho,2)),rho);
        
        vertices = [x,y];
        if sum(diff(vertices)==0) > 0
            vertices(all(diff(vertices)'==0),:) = [];
        end
        % Remove adjacent duplicates
        vertices(:,(diff(vertices(:,1))==0) & (diff(vertices(:,2))==0)) = [];
        vertices(end+1,:) = vertices(1,:);
        pt = interparc(numSamplePts,vertices(:,1),vertices(:,2),'linear');
        pt(end+1,:) = pt(1,:);
        polyshapes{i} = polyshape(pt,'Simplify',false,'KeepCollinearPoints',true);  
    else
        polyshapes{i} = nan;
    end
end



end



%------------- END CODE --------------
