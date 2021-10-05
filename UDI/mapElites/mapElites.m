function [map, percImproved, percValid, h, percImprovement] = mapElites(fitnessFunction,map,p,d)
%mapElites - Multi-dimensional Archive of Phenotypic Elites algorithm
%
% Syntax:  map = mapElites(fitnessFunction, map, p, d);
%
% Inputs:
%   fitnessFunction - funct  - returns fitness of vector of individuals
%   map             - struct - initial solutions in F-dimensional map
%   p               - struct - Hyperparameters for algorithm, visualization, and data gathering
%   d               - struct - Domain definition
%
% Outputs:
%   map    - struct - population archive
%   percImproved    - percentage of children which improved on elites
%   percValid       - percentage of children which are valid members of
%   selected classes
%   h      - [1X2]  - axes handle, data handle
%
%
% See also: createChildren, getBestPerCell, updateMap

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 02-Aug-2017

%------------- BEGIN CODE --------------

% View Initial Map
h = [];
if p.display.illu
    figure(1); clf;
    [h(1), h(2)] = viewMap(-map.fitness, d, map.edges); title('Acquisition')
    caxis([0 10]);
    figure(2); clf;
    [h(3), h(4)] = viewMap(map.confidence, d, map.edges); title('Model Variance')
    caxis([0 1]);
    figure(3); clf;
    [h(5), h(6)] = viewMap(map.confidence, d, map.edges); title('Genetic Diversity (acq)')
    caxis([0 0.5]);    
end

iGen = 1;
while (iGen < p.nGens)
    %% Create and Evaluate Children
    % Continue to remutate until enough children which satisfy geometric
    % constraints are created
    children = [];
    while size(children,1) < p.nChildren
        newChildren = createChildren(map, p, d);
        validInds = feval(d.validate,newChildren,d);
        children = [children ; newChildren(validInds,:)] ; %#ok<AGROW>
    end
    children = children(1:p.nChildren,:);
    
    % Constrain solutions to selected class(es)        
    if isfield(p, 'constraintModel') && ~isnan(p.constraints.threshold)
        isValid = applyConstraints(children, p.constraints);
        if size(isValid,1) > 1  % Multiple classes were selected
            isValid = sum(isValid);
        end
        children(~isValid,:) = []; % Remove invalid children
        percValid(iGen) = sum(isValid)/length(isValid); 
    else
        percValid(iGen) = 1;
    end
    [fitness, values] = fitnessFunction(children); %% TODO: Speed up without anonymous functions
    
    %% Add Children to Map   
    [replaced, replacement] = nicheCompete(children,fitness,map,d);
    percImproved(iGen) = length(replaced)/p.nChildren;
    map = updateMap(replaced,replacement,map,fitness,children,...
                        values,d.extraMapValues);  
    
    %% View New Map
    if p.display.illu && ~mod(iGen,p.display.illuMod)
        set(h(2),'CData',-map.fitness,...
            'AlphaData',~isnan(map.fitness))
        colormap(h(1),parula(16));
        drawnow;
        set(h(4),'CData',map.confidence,...
            'AlphaData',~isnan(map.confidence))
        colormap(h(3),parula(16));
        drawnow;
        map.diversityMap = getDiversityMap(map.genes);    
        set(h(6),'CData',map.diversityMap,...
            'AlphaData',~isnan(map.diversityMap))
        colormap(h(5),parula(16));
        drawnow;        
    end
    
iGen = iGen+1; if ~mod(iGen,100);disp([char(9) 'Illumination Generation: ' int2str(iGen)]);end;
end
if percImproved(end) > 0.05; disp('Warning: MAP-Elites finished while still making improvements ( >5% / generation )');end


%------------- END OF CODE --------------