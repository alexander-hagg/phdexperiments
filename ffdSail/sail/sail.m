function [output] = sail(p,d)
%SAIL - Surrogate Assisted Illumination Algorithm
% Main run script of SAIL algorithm
%
% Syntax:  [output] = sail(p)
%
% Inputs:
%   p  - struct - Hyperparameters for algorithm, visualization, and data gathering
%   d  - struct - Domain definition
%    * - sail with no arguments to return default parameter struct
%
% Outputs:
%    output - output struct with fields:
%               .output.model{1}             - gpModels produced by SAIL
%               .output.model{1}.trainInput  - input samples
%               .output.model{1}.trainOutput - sample results
%               .p - parameter struct (for data keeping)
%               .* - other parameters recorded throughout run
%
% Example:
%    p = sail;                                  % Load default parameters
%    p.nTotalSamples = 80;                      % Edit default parameters
%    d = af_Domain;                             % Load domain parameters
%    output = sail(d,p);                        % Run SAIL algorithm
%    predMap = createPredictionMap(output.model,p,d);
%    viewMap(predMap.fitness, p, d)             % View Prediction Map
%
% Other m-files required:
%   defaultParamSet, initialSampling, trainGP,
%   createMap, nicheCompete, updateMap, mapElites, sobolset
%
% Other submodules required: gpml-wrapper
%
% MAT-files required: if loading initial samples "d.initialSampleSource"
%
% See also: runSail, createPredictionMap, mapElites

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Nov 2016; Last revision: 03-Aug-2017

%------------- BEGIN CODE --------------
if nargin==0; output = defaultParamSet; return; end

%% 0 - Produce Initial Samples
if ~d.loadInitialSamples
    [observation, value] = initialSampling(d,p.nInitialSamples);
else
    load(d.initialSampleSource); % Contains fields: 'observation', 'value'
    randomPick = randperm(size(observation,1),p.nInitialSamples); %#ok<NODEF>
    observation = observation(randomPick,:);
    value = value(randomPick,:); %#ok<NODEF>
end
nSamples = size(observation,1);

disp(['SAIL - number of samples: ' int2str(nSamples)]);

% If a concept is defined, take the members of that concept as a seed for
% all maps
if isfield(p,'concept')
    parents = p.concept.members;
else
    parents = observation;
end

%% Acquisition Loop
trainingTime = []; trainingPredictionTime = []; illumTime = []; peTime = [];
while nSamples <= p.nTotalSamples
    %% 1 - Create Surrogate and Acquisition Function
    % Surrogate models are created from all evaluated samples, and these
    % models are used to produce an acquisition function.
    disp(['PE ' int2str(nSamples) ' | Training Acquisition Models']); tstart = tic;
    for iModel = 1:size(value,2)
        % Only retrain model parameters every 'p.trainingMod' iterations
        %if (nSamples==p.nInitialSamples || mod(nSamples,p.trainingMod*p.nAdditionalSamples));d.paramsAcq{iModel}.functionEvals = 100;else;d.paramsAcq{iModel}.functionEvals = 0;end
        d.paramsAcq{iModel}.functionEvals = 100;
        model{iModel} = feval(['train' d.paramsAcq{iModel}.name], observation, value(:,iModel), d.paramsAcq{iModel});
        %if ~mod(nSamples,p.data.mapEvalMod) % Save acquisition models for analysis later
        output.tModelsAcq{nSamples,iModel} = model{iModel};
        %end
    end
    
    % Save found model parameters and new acquisition function
    for iModel=1:size(value,2)
        if strcmp(d.paramsAcq{iModel}.name,'GP'); d.paramsAcq{iModel}.hyp = model{iModel}.hyp; end
    end
    acqFunction = feval(d.createAcqFunction, model, d);
    
    % Data Gathering
    trainingTime = [trainingTime toc(tstart)];
    
    % After final model is created no more infill is necessary
    if nSamples == p.nTotalSamples; break; end
    %% 2 - Illuminate Acquisition Map
    % A map is constructed using the evaluated samples which are evaluated
    % with the acquisition function and placed in the map as the initial
    % population. The observed samples are the seed population of the
    % 'acquisition map' which is then created by optimizing the acquisition
    % function with MAP-Elites.
    disp(['PE: ' int2str(nSamples) '| Illuminating Acquisition Map']); tstart = tic;
    % Evaluate data set with acquisition function
    try
        [fitness,predValues] = acqFunction(parents);
    catch exception
        disp(exception.identifier);
    end
    
    % Set variance normalization factor to last maximum variance and keep
    % track
    % d.prevMaxVar(end+1) = max(predValues{3}(:));
    % d.varCoef = 1./(d.prevMaxVar(end));
    
    % Place Best Samples in Map with Acquisition Fitness
    obsMap = createMap(d.featureRes, d.dof, d.extraMapValues);
    [replaced, replacement] = nicheCompete(parents, fitness, obsMap, d);
    obsMap = updateMap(replaced,replacement,obsMap,fitness,parents,...
        predValues,d.extraMapValues);
    
    % Illuminate with MAP-Elites
    [acqMap, percImproved(:,nSamples), percValid(:,nSamples)] = mapElites(acqFunction,obsMap,p,d);
    
    % Data Gathering
    illumTime = [illumTime toc(tstart)];    acqMapRecord(nSamples) = acqMap;
    %% 3 - Select Infill Samples
    % The next samples to be tested are chosen from the acquisition map: a
    % sobol sequence is used to to evenly sample the map in the feature
    % dimensions. When evaluated solutions don't converge the next bin in
    % the sobol set is chosen.
    disp(['PE: ' int2str(nSamples) '| Evaluating New Samples']); tstart = tic;
    % At first iteration initialize sobol sequence for sample selection
    if nSamples == p.nInitialSamples
        if isfield(d,'commonSobolGen')
            sobSet = d.commonSobolGen;
        else
            sobSet  = scramble(sobolset(d.nDims,'Skip',1e3),'MatousekAffineOwen');
        end
        if isfield(d,'commonSobolGenPtr')
            sobPoint = commonSobolGenPtr;
        else
            sobPoint= 1;
        end
    end
    
    newValue = nan(p.nAdditionalSamples, size(value,2)); % new values will be stored here
    noValue = any(isnan(newValue),2);
    
    while any(any(noValue))
        nNans = sum(noValue);
        nextGenes = nan(nNans,d.dof); % Create one 'blank' genome for each NAN
        
        % Identify (grab indx of NANs)
        nanIndx = 1:p.nAdditionalSamples;  nanIndx = nanIndx(noValue);
        
        % Replace with next in Sobol Sequence
        newSampleRange = sobPoint:(sobPoint+nNans-1);
        mapLinIndx = sobol2indx(sobSet,newSampleRange,d, acqMap.edges);
        
        % Replace unreachable bins
        emptyCells = isnan(acqMap.fitness(mapLinIndx));
        while any(emptyCells)
            nEmptyCells = sum(emptyCells);
            mapLinIndx(emptyCells) = sobol2indx(sobSet,sobPoint:sobPoint+nEmptyCells-1,d, acqMap.edges);
            emptyCells = isnan(acqMap.fitness(mapLinIndx));
            sobPoint = sobPoint + nEmptyCells;
        end
        
        % Pull out chosen genomes from map
        [chosenI,chosenJ] = ind2sub(d.featureRes, mapLinIndx);
        for iGenes=1:nNans
            nextGenes(iGenes,:) = acqMap.genes(chosenI(iGenes),chosenJ(iGenes),:);
        end
        
        % Precise evaluation
        if ~isempty(nextGenes)
            measuredValue = feval(d.preciseEvaluate, nextGenes, d);
            % Assign found values
            newValue(nanIndx,:) = measuredValue;
        end
        
        if p.retryInvalid
            % Check for invalid or duplicate shapes
            nanValue = any(isnan(newValue),2);
            oldDuplicate = logical(false(1,size(nanValue,1)));
            oldDuplicate(nanIndx) = any(pdist2(observation,nextGenes)==0);
            newDuplicate = logical(false(1,size(nanValue,1)));
            sampleDistances = pdist2(nextGenes,nextGenes);
            sampleDistances = sampleDistances + diag(ones(1,size(sampleDistances,1)));
            newDuplicate(nanIndx) = any(sampleDistances==0);
            noValue = nanValue | oldDuplicate' | newDuplicate';
        else
            % Do not try invalid shapes
            newValue(isnan(newValue(:,1)),:) = repmat([0 0],sum(isnan(newValue(:,1))),1);
            % We still have to skip samples from empty bins
            noValue = any(isnan(nextGenes),2);
        end
        nextObservation(nanIndx,:) = nextGenes;         %#ok<AGROW>
        sobPoint = sobPoint + length(newSampleRange);   % Increment sobol sequence for next samples
    end
    
    % Add evaluated solutions to data set
    value = cat(1,value,newValue);
    observation = cat(1,observation,nextObservation);
    nSamples  = size(observation,1);
    
    % Assign new samples to parent pool as well
    parents = cat(1,parents,nextObservation);
    
    % Data Gathering
    peTime = [peTime toc(tstart)];
end % end acquisition loop

% Create prediction map
disp(['PE ' int2str(nSamples) ' | Training Prediction Models']); tstart = tic;
for iModel = 1:size(value,2)
    d.paramsPred{iModel}.functionEvals = 100;
    modelPred{iModel} = feval(['train' d.paramsPred{iModel}.name], observation, value(:,iModel), d.paramsPred{iModel});
    output.tModelsPred{nSamples,iModel} = modelPred{iModel};
end
trainingPredictionTime = [trainingPredictionTime toc(tstart)];
[predMap(nSamples), predPercImprove(:,nSamples)] = createPredictionMap(modelPred,p,d,'featureRes',p.data.predMapRes);

% Save relevant data
d.sobPoint          = sobPoint; % Save Sobol pointer to allow experiment reuse
output.p            = p;
output.d            = d;
output.model        = model;
if exist('modelPred','var'); output.modelPred = modelPred;end
%trainingTime        = trainingTime(2:end);
output.trainingTime = trainingTime;
output.illumTime    = illumTime;
output.peTime       = peTime;
%output.time(:,4)    = trainingPredictionTime(2:end);
output.percImproved = percImproved;
if exist('predMap');output.predMap = predMap;end;
output.acqMap       = acqMapRecord;
output.unpack       = 'names = fieldnames(output); for i=1:length(names) eval( [names{i},''= output.'', names{i}, '';''] ); end';
end




