function [map, percImproved, percValid, allMaps, percFilled, stats] = illuminate(map,fitnessFunction,p,d,manifold,varargin)
%illuminate - QD with Voronoi Multi-dimensional Archive of Phenotypic Elites algorithm
%
% Syntax:  [map, percImproved, percValid, h, allMaps, percFilled] = illuminate(fitnessFunction,map,p,d,varargin)
%
% Inputs:
%   fitnessFunction - funct  - returns fitness of vector of individuals
%   map             - struct - initial solutions in F-dimensional map
%   p               - struct - Hyperparameters for algorithm, visualization, and data gathering
%   d               - struct - Domain definition
%
% Outputs:
%   map             - struct - population archive
%   percImproved    - percentage of children which improved on elites
%   percValid       - percentage of children which are valid members of selected classes
%   h               - [1X2]  - axes handle, data handle
%   allMap          - all maps created in sequence
%   percFilled      - percentage of map filled
%
%
% See also: createChildren, getBestPerCell, updateMap
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jul 2019; Last revision: 04-Jul-2019

%------------- BEGIN CODE --------------
% View Initial Map

percImproved = 0;   percValid = 0;  h = 0;  percFilled = 0;

iGen = 1;
while (iGen <= p.nGens)
    %% Create and Evaluate Children
    % Continue to remutate until enough children which satisfy geometric constraints are created
    children = [];
    while size(children,1) < p.nChildren
        newChildren = createChildren(map, p, d);
        validInds = feval(d.validate,newChildren,d);
        children = [children ; newChildren(validInds,:)] ; %#ok<AGROW>
    end
    children = children(1:p.nChildren,:);
    
    [fitness, phenotypes] = fitnessFunction(children); %% TODO: Speed up without anonymous functions
    
    %% Add Children to Map
    latent = getAELatent(manifold.model,phenotypes)';
    [replaced, replacement] = nicheCompete(children, fitness, map, d, latent);
    map = updateMap(replaced,replacement,map,fitness,children,latent,p.extraMapValues);
    
    %allMaps{iGen} = map;
    map.stats.fitnessMean(iGen) = nanmean(map.fitness);
    map.stats.fitnessTotal(iGen) = nansum(map.fitness);
    
    if ~mod(iGen,25) || iGen==1
        disp([char(9) 'Illumination Generation: ' int2str(iGen) ]);
    end
    iGen = iGen+1;
end
%if percImproved(end) > 0.05; disp('Warning: MAP-Elites finished while still making improvements ( >5% / generation )');end
end


%------------- END OF CODE --------------
