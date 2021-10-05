function showPhenotype(figHandle,d,varargin)
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

genomes = d.ranges(2)*(2*rand(1,d.dof)-1);
if nargin>2
    if ~isempty(varargin{1})
        genomes = varargin{1};
    end
end

if nargin>3
    if ~isempty(varargin{2})
        placement = varargin{2};
    end
end

selectionLabels = ones(1,size(genomes,1));
if nargin>4
    if ~isempty(varargin{3})
        selectionLabels = varargin{3};
    end
end

colors = {'red','green'};
pgon = getPhenotype(genomes,d);
for i=1:size(genomes,1)
    % Change placement if necessary
    if exist('placement','var') && ~isempty(placement)
        if isa(pgon{i},'polyshape')
            pgon{i}.Vertices = pgon{i}.Vertices + placement(i,:);
        end
    end
    if isa(pgon{i},'polyshape')
        plot(figHandle,pgon{i},'FaceColor',colors{selectionLabels(i)});
    end
    hold(figHandle,'on');
end
end

