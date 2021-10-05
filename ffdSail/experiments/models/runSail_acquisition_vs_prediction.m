%runSail - Example usage script of sail function
% Running sail without arguments will return a hyperparameter struct of
% default values. These defaults can be changed in
% /sail/defaultParamSet.m
%
% Running sail with a parameter struct as input will run the algorithm
%
%
% Other m-files required: defaultParamSet, sail, mapElites
% Other submodules required: gpml-wrapper
%
%
% See also: mapElites, sail

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander-hagg@h-brs.de
% Nov 2016; Last revision: 23-Aug-2017

%------------- BEGIN CODE --------------
clear;clc;
systemInit;
d = ffd_Domain;%d = velo_Domain; %d = af_Domain;
[d.modelPars,d.modelNames] = models_LoadCfg(d);

%%------------- EXPERIMENT LOOP --------------

p = sail;
% Adjust hyperparameters
p.nInitialSamples       = 10;
p.nAdditionalSamples    = 1;
p.nTotalSamples         = 20;
p.nChildren             = 100;
p.nGens                 = 500;

p.display.illu          = false;
p.display.illuMod       = 100;
p.display.figs          = true;

p.data.outSave          = false;
p.data.mapEval          = true;
p.data.mapEvalMod       = 10;
p.data.saveModelsMod    = 10;

% Fix initial sample set to compare models
[observation, value] = initialSampling(d,p.nInitialSamples);
d.initialSampleSource = 'initialSampleSource.mat';
save(d.initialSampleSource,'observation','value');
d.loadInitialSamples = true;

% Fix Sobol generator to compare models
d.commonSobolGen = scramble(sobolset(d.nDims,'Skip',1e3),'MatousekAffineOwen');

for modelIndex=1:length(d.modelNames)
    % % % % % % % % % % % % % % %
    disp(['Illuminating with model ' d.modelNames{modelIndex}]);
    runTime = tic;
    
    % Set parameters, including specific for HSM.
    for target=1:2;d.params{target} = feval(['params' d.modelPars{modelIndex}{1}], d.modelPars{modelIndex}{2:end});end
    
    outputs{modelIndex,1} = sail(p,d);
    disp(['Model ' d.modelNames{modelIndex} ' -- Runtime: ' seconds2human(toc(runTime))]);
    
    
    % Determine number of maps to evaluate
    iItr = p.nInitialSamples:p.data.mapEvalMod:p.nTotalSamples;
    nItr= size(iItr,2);
    
    % Get training samples and Sobol generator/counter
    observations = outputs{modelIndex,1}.model{1}.trainInput;
    values = [outputs{modelIndex,1}.model{1}.trainOutput, outputs{modelIndex,1}.model{2}.trainOutput];
    
    for secModelIndex=1:length(d.modelPars)
        output = outputs{modelIndex,1}; % For saving and reloading in ground truth script
        for i = 1:nItr            
            % Select samples for matching prediction map from SAIL run
            observation = observations(1:iItr(i),:);value = values(1:iItr(i),:);            
            disp(['Training model ' d.modelNames{secModelIndex} ' with ' int2str(size(observation,1)) ' observations']);            
            for iModel = 1:2 % lift and drag
                d.params{iModel} = feval(['params' d.modelPars{secModelIndex}{1}], d.modelPars{secModelIndex}{2:end});
                models{modelIndex,secModelIndex,i, iModel} = feval(['train' d.params{iModel}.name], observation, value(:,iModel), d.params{iModel});
            end
            disp(['Create prediction map with model ' d.modelNames{secModelIndex}]);
            output.predMap(iItr(i)) = createPredictionMap({models{modelIndex,secModelIndex,i,:}}, p, d);
        end
        output.experiment.acq = d.modelNames{modelIndex};
        output.experiment.ill = d.modelNames{secModelIndex};
        save([resultpath '/results_' d.modelNames{modelIndex} '_' d.modelNames{secModelIndex} '.mat'], 'output', 'p', 'd', '-v7.3');
    end
end





