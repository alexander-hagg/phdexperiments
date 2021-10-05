function [predMap,fitModelPred,featModelsPred,allMaps,trueFilled] = sphen(sampleSet,p,d,varargin)
%SPHEN - Surrogate-Assisted Phenotypic Niching
% Main run script of SPHEN algorithm
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Apr 2020; Last revision: 09-Apr-2020

axHandle = [];
if p.display.illu
    if nargin > 3; figID = varargin{1};else; figID = 1; end
    f = figure(figID);clf(f);
    axHandle = gca;
    figIllumination(1) = figure(figID+1); fitHandle = gca;
    figIllumination(2) = figure(figID+2); featHandle = gca;
    figIllumination(3) = figure(figID+3); trueFillHandle = gca;
    figIllumination(4) = figure(figID+4); predMapHandle = gca;
end

observation = sampleSet.samples;
trueFitness = sampleSet.fitness;
trueFeatures = sampleSet.features;

p.numInitSamples = size(observation,1);
nSamples = p.numInitSamples;

p.infill.modelParams = paramsGP(d.dof);
predMapCtr = 1;
allMaps{predMapCtr} = [];
allACQMaps{predMapCtr} = [];
trueFilled = 0;

while nSamples <= p.infill.nTotalSamples
    %% 1 - Create Surrogate and Acquisition Function
    % Surrogate models are created from all evaluated samples, and these
    % models are used to produce an acquisition function. 
    disp(['SPHEN: PE: ' int2str(nSamples) ' | Training Models']);
    [fitModel,featModels] = trainModels(observation,trueFitness,trueFeatures,p.infill.modelParams);
    acqFunction = @(x) ucb(x, [fitModel, featModels{:}], d, p);
    
    % After final model is created no more infill is necessary
    if nSamples >= p.infill.nTotalSamples; break; end
    
    %% 2 - Illuminate Acquisition Map
    % A map is constructed using the samples which are evaluated
    % with the acquisition function and placed in the map as the initial
    % population. The observed samples are the seed population of the
    % 'acquisition map' which is then created by optimizing the acquisition
    % function with QD.
    
    if p.infill.injectTrueValues
        % Inject true values of samples into acquisition map
        fitness = trueFitness;
        features = trueFeatures;
    else
        % Evaluate data set with acquisition function
        try
            [fitness,features] = acqFunction(observation);
        catch exception
            disp(exception.identifier);
        end
    end
    
    % Place Best Samples in Map with Acquisition Fitness
    acqMap                                              = createMap(d, p);
    [replaced, replacement]                             = nicheCompete(observation, fitness, acqMap, d, p, features);
    acqMap                                              = updateMap(replaced,replacement,acqMap,fitness,observation,features);
    
    % Illuminate with QD
    disp([char(9) 'Creating Acquisition Map']);
    acqMap = illuminate(acqMap,acqFunction,p,d,axHandle);
        
    %% 3 - Select Infill Samples
    % The next samples to be tested are chosen from the acquisition map: a
    % sobol sequence is used to to evenly sample the map in the feature
    % dimensions. When evaluated solutions don't converge the next bin in
    % the sobol set is chosen.
    disp([char(9) 'Evaluating New Samples']);
    
    % At first iteration initialize sobol sequence for sample selection
    if nSamples == p.numInitSamples
        sobSet  = scramble(sobolset(numel(acqMap.resolution),'Skip',1e3),'MatousekAffineOwen');
        sobPoint= 1;
    end
    
    % Choose new samples and evaluate them for new observations
    nMissing = p.infill.nAdditionalSamples; newValue = []; newSample = []; newFeatures = [];
    while nMissing > 0
        % Evenly sample solutions from acquisition map
        newSampleRange      = sobPoint:(sobPoint+100*p.infill.nAdditionalSamples)-1;
        [~,binIndx]        = sobol2indx(sobSet,newSampleRange, size(acqMap.fitness), acqMap.edges);
        for iGenes=1:size(binIndx,1)
            indPool(iGenes,:) = acqMap.genes(binIndx(iGenes,1),binIndx(iGenes,2),:); % TODO: more than 2 dims
        end
        % Remove repeats and nans (empty bins)
        indPool = setdiff(indPool,observation,'rows','stable');
        indPool = indPool(~any(isnan(indPool),2),:);
        
        % Evaluate enough of these valid solutions to get your initial sample set
        peFunction = @(x) feval(d.fitfun,x);
        [foundSample, foundValue, foundFeatures, nMissing] = getValidInds(indPool, peFunction, nMissing, p.infill.retryInvalid);
        newSample = [newSample; foundSample]; %#ok<*AGROW>
        newValue  = [newValue ; foundValue ]; %#ok<*AGROW>
        newFeatures = [newFeatures ; foundFeatures ]; %#ok<*AGROW>
        
        % Advance sobol sequence
        sobPoint = sobPoint + p.infill.nAdditionalSamples + 1;
    end    
    
    % Add evaluated solutions to data set
    trueFitness = cat(1,trueFitness,newValue);
    observation = cat(1,observation,newSample);
    trueFeatures = cat(1,trueFeatures,newFeatures);
    nSamples  = size(observation,1);
    
    
    
    if p.infill.display.illu% && (~mod(iGen,p.display.illuMod) || (iGen==p.nGens))
        trueFilled = visualizeStats(acqMap,d,nSamples,sampleSet,p,trueFeatures,newFeatures,trueFilled,fitHandle,featHandle,trueFillHandle);
    end
    
    if p.infill.intermediateMaps && ~mod(predMapCtr,p.infill.intermediateMapSteps)
        disp([char(9) 'Creating Prediction Map']);
        allMaps{predMapCtr} = createPredictionMap([fitModel, featModels{:}],p,d);        
        allACQMaps{predMapCtr} = acqMap;
        save('allMapsSPHEN.mat','allMaps','allACQMaps','fitModel','featModels','trueFitness','observation','trueFeatures');
        predMapCtr = predMapCtr + 1;
    end
    
end % end acquisition loop


% Create prediction map
disp([char(9) 'Creating Final Prediction Map']);
[fitModelPred,featModelsPred] = trainModels(observation,trueFitness,trueFeatures,p.infill.modelParams);
[predMap] = createPredictionMap([fitModelPred, featModelsPred{:}],p,d);

if p.display.illu
    viewMap(predMap,d,'fitness',predMapHandle);caxis([0 1]);
end

end

function trueFilled = visualizeStats(map,d,nSamples,sampleSet,p,trueFeatures,newFeatures,trueFilled,fitHandle,featHandle,trueFillHandle)

tGenes = reshape(map.genes,[],16); 
nanVals = any(isnan(tGenes'));
tGenes(nanVals,:) = [];
[tFitness,tFeatures] = d.fitfun(tGenes);

tFeatures(:,1) = (tFeatures(:,1)-d.featureMin(1))./(d.featureMax(1)-d.featureMin(1));
tFeatures(:,2) = (tFeatures(:,2)-d.featureMin(2))./(d.featureMax(2)-d.featureMin(2));
tFeatures(tFeatures>1) = 1; tFeatures(tFeatures<0) = 0;

truemap = createMap(d, p);
[replaced, replacement] = nicheCompete(tGenes, tFitness, truemap, d, p, tFeatures);
truemap = updateMap(replaced,replacement,truemap,tFitness,tGenes,tFeatures);

trueFilled(end+1) = 100*sum(~isnan(truemap.fitness(:)))/numel(truemap.fitness(:));

viewMap(truemap,d,'fitness',fitHandle);
caxis(fitHandle,[0 1]);
title(fitHandle,['True - ' int2str(nSamples) ' samples - ' num2str(trueFilled(end)) '%fill']);

hold(featHandle,'off');
scatter(featHandle,sampleSet.features(:,1),sampleSet.features(:,2),32,'r','filled');
hold(featHandle,'on');
scatter(featHandle,trueFeatures(p.numInitSamples+1:end-size(newFeatures,1),1),trueFeatures(p.numInitSamples+1:end-size(newFeatures,1),2),32,'k','filled');
scatter(featHandle,newFeatures(:,1),newFeatures(:,2),16,'b','filled');
axis(featHandle,[d.featureMin(1) d.featureMax(1) d.featureMin(2) d.featureMax(2)]);
axis(featHandle,'square');

hold(trueFillHandle,'off');
plot(trueFillHandle,trueFilled);
title(trueFillHandle,'True Filled %');
trueFillHandle.YAxis.Limits = [0 100];
grid(trueFillHandle,'on');
trueFillHandle.XTick = 1:numel(trueFilled);
step = (nSamples-numel(sampleSet.fitness))/(numel(trueFilled)-1);
trueFillHandle.XTickLabel = numel(sampleSet.fitness):step:nSamples;
xlabel(trueFillHandle,'# Samples');
drawnow;

end
