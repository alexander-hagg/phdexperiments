function [map, percImproved, percValid, allMaps, percFilled] = illuminate(constraints,p,d,varargin)
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
simspaceaxes = []; if nargin > 3; simspaceaxes = varargin{1};end
if nargin > 4; figHandleMap = varargin{2};else;f=figure(1);clf(f);figHandleMap = axes; end
if nargin > 5; figHandleTotalFit = varargin{3};else;f=figure(2);clf(f);figHandleTotalFit = axes;end
if nargin > 6; figHandleMeanDrift = varargin{4};else;f=figure(3);clf(f);figHandleMeanDrift = axes;end

fitnessFunction = @(x) objective(x,d.fitfun,[],p.penaltyWeight,p.driftThreshold,simspaceaxes);
if ~isempty(constraints)
    disp(['Adding constraints to the fitness function']);
    fitnessFunction = @(x) objective(x,d.fitfun,constraints,p.penaltyWeight,p.driftThreshold,simspaceaxes);
    initSamples = [];
    for it1=1:length(constraints)
        initSamples = [initSamples; constraints{it1}.members];
    end
else
    disp(['Initializing space filling sample set']);
    sobSequence         = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');
    sobPoint            = 1;
    initSamples         = range(d.ranges).*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(1);
end
[fitness, values, phenotypes]       = fitnessFunction(initSamples); %

map                                 = createMap(d, p);
[replaced, replacement, features]   = nicheCompete(initSamples, fitness, phenotypes, map, d, p);
map                                 = updateMap(replaced,replacement,map,fitness,initSamples,values,features,p.extraMapValues);

if p.display.illu; viewMap(map,d,figHandleMap); title(figHandleMap,'Fitness'); caxis(figHandleMap,[0 1]);drawnow;end

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
title(figHandleTotalFit,['Total Fitness (QD Fitness)']);

if sum(driftMean(:)) > 0
    cla(figHandleMeanDrift);
    plot(figHandleMeanDrift,driftMean,'LineWidth',2);
    axis(figHandleMeanDrift,[0 iGen 0 1]);
    grid(figHandleMeanDrift,'on');
    title(figHandleMeanDrift,['User Selection Drift']);    
end

drawnow;
end
%------------- END OF CODE --------------
