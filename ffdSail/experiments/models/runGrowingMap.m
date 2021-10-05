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
clear;clc;cd('../..');systemInit;
d = ffd_Domain; %d = velo_Domain; %d = af_Domain;
[d.modelParsAcq, d.modelParsPre] = models_LoadCfg(d);

%%------------- EXPERIMENT LOOP --------------

p = sail;
% Adjust hyperparameters
p.nInitialSamples       = 10;
p.nAdditionalSamples    = 10;
p.nTotalSamples         = 50;
p.nChildren             = 100;
p.nGens                 = 500;
p.retryInvalid          = true;

p.display.illu          = false;
p.display.illuMod       = 100;
p.display.figs          = true;

p.data.outSave          = false;
p.data.mapEval          = true;
p.data.mapEvalMod       = 50;
p.data.mapEvalSteps     = [p.nInitialSamples p.data.mapEvalMod:p.data.mapEvalMod:p.nTotalSamples]; % Default:

% Fix initial sample set and Sobol generator for fair model comparison
[observation, value] = initialSampling(d,p.nInitialSamples);
d.initialSampleSource = 'initialSampleSource.mat'; save(d.initialSampleSource,'observation','value'); d.loadInitialSamples = true;
d.commonSobolGen = scramble(sobolset(d.nDims,'Skip',1e3),'MatousekAffineOwen');

for modelIndex=1:length(d.modelParsAcq)
    % % % % % % % % % % % % % % %
    disp(['Illuminating with model ' d.modelParsAcq{modelIndex}{1}]);
    runTime = tic;
    
    % Set acquisition and prediction model parameters
    for target=1:2
        d.params{target}        = feval(['params' d.modelParsAcq{modelIndex}{2}], d.modelParsAcq{modelIndex}{3:end});
        d.paramsPred{target}    = feval(['params' d.modelParsPre{modelIndex}{2}], d.modelParsPre{modelIndex}{3:end});
    end
    
    output = sail(p,d);
    
    disp(['Acquisition ' d.modelParsAcq{modelIndex}{1} ' , Prediction ' d.modelParsPre{modelIndex}{1} '-- Runtime: ' seconds2human(toc(runTime))]);
    save([resultpath '/results_' d.modelParsAcq{modelIndex}{1} '_' d.modelParsPre{modelIndex}{1} '.mat'], 'output', 'p', 'd', '-v7.3');
end





