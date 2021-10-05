function map = createMap(d, p, varargin)
%createMap - Defines map struct and feature space cell divisions
%
% Syntax:  map = createMap(d, p, varargin)
%
% Inputs:
%    d - struct - Domain description struct
%    p - struct - QD configuration struct
%
% Optional Inputs:
%    varargin - extra map values
%
% Outputs:
%    map  - struct with [M(1) X M(2)...X M(N)] matrices for fitness, etc
%       edges               - {1XN} cell of partitions for each dimension
%       fitness, drag, etc  - [M(1) X M(2) X M(N)] matrices of scalars
%       genes               - [M(1) X M(2) X M(N) X genomeLength]
%
% Example: 
%   map = createMap([10 5], 3); % 10 X 5 map of genomes with 3 parameters
%   OR
%   extraMapValues = {'cD','cL'};
%   map = createMap(d.featureRes, d.dof, extraMapValues)
%
% Author: Adam Gaier, Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de, alexander.hagg@h-brs.de
% Jun 2016; Last revision: 15-Aug-2019

%------------- BEGIN CODE --------------

for i=1:length(p.featureResolution)
    edges{i} = linspace(0,1,p.featureResolution(i)+1); %#ok<AGROW>
end
map.edges = edges;

blankMap     = NaN(p.featureResolution,'double');
map.fitness  = blankMap;
map.genes       = repmat(blankMap,[1 1 d.dof]); %#ok<REPMAT>

% Phenotypic Features
map.features = repmat(blankMap,[1 1 2]); %#ok<REPMAT>

if ~isempty(varargin)
    for iValues = 1:length(varargin{1})
        eval(['map.' varargin{1}{iValues} ' = blankMap;']);
    end
end


%------------- END OF CODE --------------