function [flatbooleanMap,polyshapes] = getPhenotypeFFD(genomes,base,varargin)
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
tol = 0.01;
resolution = 64; if nargin>2;resolution = varargin{1};end
nGenomes = 1; if size(genomes,2)>1; nGenomes = size(genomes,1); end

for i=1:nGenomes
    if size(genomes,2)>1; genome = genomes(i,:); else;genome = genomes';end
    if ~isnan(genome(1))
        [theta,rho] = cart2pol(base(1:end/2)',base(end/2+1:end)');
        rho         = rho .* (genome(1:end/2)');
        theta       = theta + genome(end/2+1:end)';
        [x,y]       = pol2cart(repmat(theta,1,size(rho,2)),rho);
        
        polyshapes{i} = polyshape(x,y,'Simplify',false,'KeepCollinearPoints',true);
        pixelCoordinates{i} = ceil(polyshapes{i}.Vertices*(resolution)/2)+resolution/2;
        pixelCoordinates{i}(all(isnan(pixelCoordinates{i})'),:) = [];
        booleanMap{i} = poly2mask(pixelCoordinates{i}(:,1),pixelCoordinates{i}(:,2),resolution,resolution);
    else
        polyshapes{i} = nan;
    end
end

flatbooleanMap = cat(3,[],booleanMap{:});
flatbooleanMap = reshape(flatbooleanMap,[],size(flatbooleanMap,3));
flatbooleanMap = flatbooleanMap';

end



%------------- END CODE --------------