function map = updateMap(replaced,replacement,map,fitness,genes,features,extraMapValues,values)
%updateMap - Replaces all values in a set of map cells
%
% Syntax:  map = updateMap(replaced,replacement,map,fitness,genes,values,features,extraMapValues)
%
% Inputs:
%   replaced    - [1XM]  - linear index of map cells to be replaced
%   replacement - [1XM]  - linear index of children values to place in map
%   map         - struct - population archive
%   fitness     - [1XN]  - Child fitness
%   genes       - [NXD]  - Child genomes
%   values      - [1XN]  - extra values of interest, e.g. 'cD'
%   features    - [NxnumFeatures]
%   extraMapValues
%
% Outputs:
%   map         - struct - population archive
%
%
% See also: createMap, nicheCompete

% Author: Adam Gaier, Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de, alexander.hagg@h-brs.de
% Jun 2016; Last revision: 15-Aug-2019

%------------- BEGIN CODE --------------

% Assign Fitness
map.fitness (replaced) = fitness (replacement)';

% Assign Genomes and Features
[r,c] = size(map.fitness);
[replacedI,replacedJ] = ind2sub([r c], replaced);
for iReplace = 1:length(replaced)
    map.genes(replacedI(iReplace),replacedJ(iReplace),:) = ...
        genes(replacement(iReplace),:) ;         
    map.features(replacedI(iReplace),replacedJ(iReplace),:) = ...
        features(replacement(iReplace),:) ;             
end

% Assign Miscellaneous Map values
if exist('extraMapValues','var') && ~isempty(extraMapValues)
    for iValues = 1:length(extraMapValues)
        eval(['map.' extraMapValues{iValues} '(replaced) = values{' int2str(iValues) '}(replacement);']);        
    end
end

%------------- END OF CODE --------------