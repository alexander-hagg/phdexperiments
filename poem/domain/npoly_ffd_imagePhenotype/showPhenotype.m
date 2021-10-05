function figHandle = showPhenotype(genomes,d,varargin)
%showPhenotype - Either show an example phenotype, or, when given, show
%                multiple phenotypes that are positioned on predefined placement positions.
%                Yes, this visualization script does too many things at the same time.
%
% Syntax:   showPhenotype(figHandle,d,varargin)
%
% Inputs:
%    figHandle      - [1] - figure handle
%    d              - struct - Domain description struct
%
% Optional Inputs:
%    varargin{1}    - [NxM] - genomes
%    varargin{2}    - [Nx2] - placement
%    varargin{3}    - [Nx1] - selection labels (1 or 2)
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Aug 2019; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------
nShapes = size(genomes,1);
xPos = 0:ceil(sqrt(nShapes))-1; [X,Y] = ndgrid(xPos,xPos);
placement = [X(:) Y(:)]; placement = placement(1:nShapes,:); placement = 2*placement;
if nargin>2
    figHandle = varargin{1};
else
    figHandle = figure;
end
if nargin>3
    if ~isempty(varargin{2})
        placement = varargin{2};
    end
end

if nargin>4
    clrs = varargin{3};
end

pgon = d.getPhenotype(genomes);
figure(figHandle);hold off;
for i=1:nShapes
    % Change placement if necessary
    if exist('placement','var') && ~isempty(placement)
        pgon{i}.Vertices = pgon{i}.Vertices + placement(i,:);
    end
    if isa(pgon{i},'polyshape')
        plot(pgon{i},'FaceColor',clrs(i,:));
    end
    hold('on');
end

axis([-5 max(placement(:))+5 -5 max(placement(:))+5]);
end

