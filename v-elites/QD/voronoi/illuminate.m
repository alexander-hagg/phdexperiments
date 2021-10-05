function [map, percImproved, percValid, allMaps, percFilled] = illuminate(fitnessFunction,map,p,d,varargin)
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
if nargin > 4; figHandleMap = varargin{1};else;figHandleMap = figure(1);end
if nargin > 5; figHandleTotalFit = varargin{2};else;figHandleTotalFit = figure(2);end
if nargin > 6; figHandleMeanDrift = varargin{3};else;figHandleMeanDrift = figure(3);end

if p.display.illu
    cla(figHandleMap);cla(figHandleTotalFit);cla(figHandleMeanDrift);
    viewMap(map,d,figHandleMap); title(figHandleMap,'Fitness'); caxis(figHandleMap,[0 1]);
    drawnow;
end

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
    
    [fitness, values, phenotypes] = fitnessFunction(children); %% TODO: Speed up without anonymous functions
    
    %% Add Children to Map
    [replaced, replacement, features] = nicheCompete(children, fitness, phenotypes, map, d);
    map = updateMap(replaced,replacement,map,fitness,children,values,features,p.extraMapValues);
    %    percImproved(iGen) = length(replaced)/p.nChildren;
    
    allMaps{iGen} = map;
    %    percFilled(iGen) = sum(~isnan(map.fitness(:)))/(size(map.fitness,1)*size(map.fitness,2));
    fitnessMean(iGen) = nanmean(map.fitness);
    fitnessTotal(iGen) = nansum(map.fitness);
    driftMean(iGen) = nanmean(map.drift);
    
    %% View New Map
    if p.display.illu && (~mod(iGen,p.display.illuMod) || (iGen==p.nGens))
        visualizeStats(figHandleMap,figHandleTotalFit,figHandleMeanDrift,iGen,p.nGens,numel(map.fitness),fitnessTotal,driftMean,map,d)
    end
    if ~mod(iGen,25) || iGen==1
        disp([char(9) 'Illumination Generation: ' int2str(iGen) ]);
    end
    iGen = iGen+1;
end
%if percImproved(end) > 0.05; disp('Warning: MAP-Elites finished while still making improvements ( >5% / generation )');end
end

function visualizeStats (figHandle,figHandleTotalFit,figHandleMeanDrift,iGen,nGens,numElements,fitnessTotal,driftMean,map,d)
cla(figHandle);
viewMap(map,d,figHandle);
title(figHandle,['Fitness Gen ' int2str(iGen) '/' int2str(nGens)]);
caxis(figHandle,[0 1]);

cla(figHandleTotalFit);
plot(figHandleTotalFit,fitnessTotal./numElements,'LineWidth',2);
axis(figHandleTotalFit,[0 iGen 0 1]);
grid(figHandleTotalFit,'on');

if sum(driftMean(:)) > 0
    cla(figHandleMeanDrift);
    plot(figHandleMeanDrift,driftMean,'LineWidth',2);
    axis(figHandleMeanDrift,[0 iGen 0 1]);
    grid(figHandleMeanDrift,'on');
end

drawnow;
end
%------------- END OF CODE --------------
