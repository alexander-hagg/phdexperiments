function [predMap,model,allMaps] = sail(sampleSet,p,d,varargin)
%SAIL - Surrogate Assisted Illumination Algorithm
% Main run script of SAIL algorithm
%
% Author: Adam Gaier, Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de, alexander.hagg@h-brs.de
% Nov 2016; Last revision: 09-Apr-2020

if p.display.illu
    if nargin > 6; figHandleAcqMap = varargin{4};else;f=figure(4);clf(f);figHandleAcqMap = axes; end
    title(figHandleAcqMap,'Acquisition Fcn'); drawnow;
end

observation = sampleSet.samples;
trueFitness = sampleSet.fitness;

p.numInitSamples = size(observation,1);
nSamples = p.numInitSamples;

p.infill.modelParams = paramsGP(d.dof);
predMapCtr = 1;

while nSamples <= p.infill.nTotalSamples
    %% 1 - Create Surrogate and Acquisition Function
    % Surrogate models are created from all evaluated samples, and these
    % models are used to produce an acquisition function.
    disp(['SAIL: PE: ' int2str(nSamples) ' | Training Model']);
    model = trainGP(observation,trueFitness,p.infill.modelParams);
    acqFunction = @(x) ucb(x, model, d, p);
    
    % After final model is created no more infill is necessary
    if nSamples >= p.infill.nTotalSamples; break; end
    
    %% 2 - Illuminate Acquisition Map
    % A map is constructed using the samples which are evaluated
    % with the acquisition function and placed in the map as the initial
    % population. The observed samples are the seed population of the
    % 'acquisition map' which is then created by optimizing the acquisition
    % function with QD.
    
    % Evaluate data set with acquisition function
    try
        [fitness,features] = acqFunction(observation);
    catch exception
        disp(exception.identifier);
    end
    
    
    % Place Best Samples in Map with Acquisition Fitness
    acqMap                                              = createMap(d, p);
    [replaced, replacement, features]                   = nicheCompete(observation, fitness, acqMap, d, p, features);
    acqMap                                              = updateMap(replaced,replacement,acqMap,fitness,observation,features);
        
    % Illuminate with QD
    disp([char(9) 'Creating Acquisition Map']);
    acqMap = illuminate(acqMap,acqFunction,p,d);
    
    %% 3 - Select Infill Samples
    % The next samples to be tested are chosen from the acquisition map: a
    % sobol sequence is used to to evenly sample the map in the feature
    % dimensions. When evaluated solutions don't converge the next bin in
    % the sobol set is chosen.
    disp([char(9) ' Evaluating New Samples']);
    
    % At first iteration initialize sobol sequence for sample selection
    if nSamples == p.numInitSamples
        if isfield(d,'commonSobolGen')
            sobSet = d.commonSobolGen;
        else
            sobSet  = scramble(sobolset(numel(acqMap.resolution),'Skip',1e3),'MatousekAffineOwen');
        end
        if isfield(d,'commonSobolGenPtr'); sobPoint = commonSobolGenPtr; else; sobPoint= 1; end
    end
    
    newValue = nan(p.infill.nAdditionalSamples, size(fitness,2)); % new values will be stored here
    noValue = any(isnan(newValue),2);
    
    while any(any(noValue))
        nNans = sum(noValue);
        nextGenes = nan(nNans,d.dof); % Create one 'blank' genome for each NAN
        
        % Identify (grab indx of NANs)
        nanIndx = 1:p.infill.nAdditionalSamples;  nanIndx = nanIndx(noValue);
        
        % Replace with next in Sobol Sequence
        newSampleRange = sobPoint:(sobPoint+nNans-1);
        mapLinIndx = sobol2indx(sobSet, newSampleRange, acqMap.resolution, acqMap.edges);
        
        % Replace unreachable bins
        emptyCells = isnan(acqMap.fitness(mapLinIndx));
        while any(emptyCells)
            nEmptyCells = sum(emptyCells);
            mapLinIndx(emptyCells) = sobol2indx(sobSet,sobPoint:sobPoint+nEmptyCells-1, acqMap.resolution, acqMap.edges);
            emptyCells = isnan(acqMap.fitness(mapLinIndx));
            sobPoint = sobPoint + nEmptyCells;
        end
        
        % Pull out chosen genomes from map
        [chosenI,chosenJ] = ind2sub(acqMap.resolution, mapLinIndx);
        for iGenes=1:nNans
            nextGenes(iGenes,:) = acqMap.genes(chosenI(iGenes),chosenJ(iGenes),:);
        end
        
        % Precise evaluation
        if ~isempty(nextGenes)
            measuredValue = feval(d.fitfun,nextGenes);
            % Assign found values
            newValue(nanIndx,:) = measuredValue;
        end
        
        % List invalid shapes
        nanValue = any(isnan(newValue),2);
        if p.infill.retryInvalid
            % Check for duplicate shapes      
            oldDuplicate = logical(false(1,size(nanValue,1)));
            oldDuplicate(nanIndx) = any(pdist2(observation,nextGenes)==0);
            newDuplicate = logical(false(1,size(nanValue,1)));
            sampleDistances = pdist2(nextGenes,nextGenes) + diag(ones(1,size(nextGenes,1)));
            newDuplicate(nanIndx) = any(sampleDistances==0);
            noValue = nanValue | oldDuplicate' | newDuplicate';
        else % Do not retry invalid shapes            
            newValue(nanValue,:) = repmat([5],sum(nanValue),1);
            % Skip samples from empty bins
            noValue = any(isnan(nextGenes),2);
        end
        nextObservation(nanIndx,:) = nextGenes;         %#ok<AGROW>
        sobPoint = sobPoint + length(newSampleRange);   % Increment sobol sequence for next samples
    end
    
    % Add evaluated solutions to data set
    trueFitness = cat(1,trueFitness,newValue);
    observation = cat(1,observation,nextObservation);
    nSamples  = size(observation,1);
    
    if p.infill.display.illu% && (~mod(iGen,p.display.illuMod) || (iGen==p.nGens))
        visualizeStats(figHandleAcqMap,acqMap,d,nSamples);
    end
    
    if p.infill.intermediateMaps
        disp([char(9) 'Training Prediction Models']);
        allMaps{predMapCtr} = createPredictionMap(model,p,d);
        predMapCtr = predMapCtr + 1;
    end
    
end % end acquisition loop


% Create prediction map
disp([char(9) 'Creating Final Prediction Map']);
[predMap] = createPredictionMap(model,p,d);

if p.display.illu
    viewMap(predMap,d);
end

end

function visualizeStats(figHandle,map,d,nSamples)
cla(figHandle);
caxis(figHandle,[0 max(map.fitness(:))]);
viewMap(map,d,'fitness',figHandle);
title(figHandle,['Acquisition Fcn w. ' int2str(nSamples) ' samples']);
drawnow;
end
