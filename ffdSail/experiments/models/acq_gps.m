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
domainname = 'PARSEC';
cd('../..');systemInit;

%%
if strcmp(domainname,'PARSEC');d = parsec_Domain(tmpdir); end
if strcmp(domainname,'FOILFFD');d = ffd_Domain(tmpdir); end
d.fitnessPenalty_Area = false;
[d.modelParsAcq, d.modelParsPre] = models_LoadCfg(d);

%%------------- EXPERIMENT LOOP --------------
p = sail;
% Adjust hyperparameters
p.nInitialSamples       = 50;
p.nAdditionalSamples    = 10;
p.nTotalSamples         = 500;
p.nChildren             = 100;
p.nGens                 = 500;
p.retryInvalid          = true;

p.display.illu          = false;
p.display.illuMod       = 100;
p.display.figs          = true;

p.data.outSave          = false;
p.data.mapEval          = true;
p.data.mapEvalMod       = 50;
p.data.mapEvalSteps     = [p.data.mapEvalMod:p.data.mapEvalMod:p.nTotalSamples];

% Fix initial sample set to compare models
d.loadInitialSamples = false;
if d.loadInitialSamples
    [observation, value] = initialSampling(d,p.nInitialSamples);
    d.initialSampleSource = 'initialSampleSource.mat';
    save(d.initialSampleSource,'observation','value');
    
    % Fix Sobol generator to compare models
    d.commonSobolGen = scramble(sobolset(d.nDims,'Skip',1e3),'MatousekAffineOwen');
end



% % % % % % % % % % % % % % %
disp(['Illuminating with model ' d.modelParsAcq{1}]);    
runTime = tic;
    
% Set acquisition and prediction model parameters
for target=1:2
    d.params{target}        = feval(['params' d.modelParsAcq{2}], d.modelParsAcq{3:end});
    d.paramsPred{target}    = feval(['params' d.modelParsPre{2}], d.modelParsPre{3:end});
end
    
output = sail(p,d);
    
disp(['Acquisition ' d.modelParsAcq{1} ' , Prediction ' d.modelParsPre{1} '-- Runtime: ' seconds2human(toc(runTime))]);
save([resultpath '/results_' d.modelParsAcq{1} '_' d.modelParsPre{1} '.mat'], 'output', 'p', 'd', '-v7.3');






