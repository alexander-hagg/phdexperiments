function [map, percImproved, percValid, allMaps, percFilled, fitnessMean, fitnessTotal] = illuminate(map,fitnessFunction,p,d,varargin)
%illuminate - QD with Multi-dimensional Archive of Phenotypic Elites algorithm
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
% Author: Adam Gaier, Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de, alexander.hagg@h-brs.de
% Jun 2016; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------
if p.display.illu
    if nargin > 4; axHandle = varargin{1};else;f=figure(99);clf(f);axHandle = gca; end
    viewMap(map,d,'fitness',axHandle); title(axHandle,'Fitness'); caxis(axHandle,[0 1]);drawnow;        
end

percFilledInit = sum(~isnan(map.fitness(:)))/(size(map.fitness,1)*size(map.fitness,2));
disp([char(9) 'Initial Map Coverage: ' num2str(100*percFilledInit) '%']);
    
iGen = 1;
while (iGen <= p.nGens)
    %% Create and Evaluate Children
    % Continue to remutate until enough children which satisfy geometric constraints are created
    children = []; parentMapIDs = [];
    while size(children,1) < p.nChildren
        newChildren = createChildren(map, p, d);
        validInds = feval(d.validate,newChildren,d);
        children = [children ; newChildren(validInds,:)] ; %#ok<AGROW>
    end
    children = children(1:p.nChildren,:);
    
    percValid(iGen) = 1;
    [fitness, features] = fitnessFunction(children);
    
    %% Add Children to Map
    [replaced, replacement, features] = nicheCompete(children,fitness,map,d,p,features);
    percImproved(iGen) = length(replaced)/p.nChildren;
    map = updateMap(replaced,replacement,map,fitness',children,features);    
    
    allMaps{iGen} = map;
    percFilled(iGen) = sum(~isnan(map.fitness(:)))/(size(map.fitness,1)*size(map.fitness,2));
    fitnessMean(iGen) = nanmean(map.fitness(:));
    fitnessTotal(iGen) = nansum(map.fitness(:));
    
    %% View New Map
    if p.display.illu && (~mod(iGen,p.display.illuMod) || (iGen==p.nGens))
        viewMap(map,d,'fitness',axHandle);
        title(axHandle,['Fitness Gen ' int2str(iGen) '/' int2str(p.nGens)]);
        caxis(axHandle,[0 1]);
        drawnow;        
    end
    if ~mod(iGen,256) || iGen == p.nGens
        disp([char(9) 'Illumination Generation: ' int2str(iGen) ' - Map Coverage: ' num2str(100*percFilled(iGen)) '% - Improvement: ' num2str(100*percImproved(iGen))]);
    end
    iGen = iGen+1;
end

% Warning if MAP-Elites hasn't converged yet (limit set in defaultParamSet)
if percImproved(end) > p.convergeLimit; disp(['Warning: MAP-Elites finished while still making improvements ( >' num2str(p.convergeLimit*100) '% / generation )']);end
end


%------------- END OF CODE --------------
