%RUNPRODIGI - Example usage script of PRODIGI: Prototype Discovery
%Generative Iterations
%
% Other submodules required: sail, gpml
%

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2018; Last revision: 25-Jan-2018

%------------- BEGIN CODE --------------
% Clean up workspace and add relevant files to path
clear;

% Domain
domainname = 'MIRROR'; % domainname = 'FOILFFD';
systemInit;
%d = mirror_Domain('hpc', false, 'lowres', true); %d = ffd_Domain('/scratch/ahagg2s/sailCFD/tmp');
d = mirror_Domain('hpc', true, 'lowres', true, 'nCases', 2, 'runFolder', '/home/ahagg2s/prodigi');

% Algorithm hyperparameters
p = produqd;                                                % load default hyperparameters

% ADJUSTMENT FOR TESTING
p.qd.nGens              = 512;
p.qd.nInitialSamples    = 20;
p.qd.nAdditionalSamples = 10;
p.qd.nTotalSamples      = 40;
p.numIterations         = 2;

% Edit hyperparameters
p.qd.data.mapEval           = true;                         % produce intermediate prediction maps
p.qd.data.mapEvalMod        = p.qd.nAdditionalSamples;      % every n samples

% Set acquisition and prediction model parameters
for target=1:size(d.modelParsAcq,1)
    d.paramsAcq{target}     = feval(['params' d.modelParsAcq{target,2}], d.modelParsAcq{target,3:end});
    d.paramsPred{target}    = feval(['params' d.modelParsPre{target,2}], d.modelParsPre{target,3:end});
end
p.qd.data.predMapRes    = d.featureRes;

%% Run PRODIGI
runTime = tic;
output = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);

%% Get Ground Truth
if strcmp(domainname,'FOILFFD')
    for it=1:length(output)
        predMap = getGroundTruthMapArray(output{it}.sail.predMap(end), d, 'dummy', false);
        output{it}.sail.predMap(end).fitness_true = predMap.fitness_true;
        output{it}.sail.predMap(end).cD_true = predMap.cD_true;
        output{it}.sail.predMap(end).cL_true = predMap.cL_true;
    end
elseif strcmp(domainname,'MIRROR')
    for it=1:length(output)
        predMap = getGroundTruthMapArrayMirror(output{it}.sail.predMap(end), d, 'dummy', false);
        output{it}.sail.predMap(end).fitness_true = predMap.fitness_true;
        output{it}.sail.predMap(end).dragForce_true = predMap.dragForce_true;
    end
end
%% Visualization
%fig(1) = figure(1);
%[~,~,cH] = viewMap(output{end}.sail.predMap(end).fitness,d);caxis([0.4 1]);cH.Label.String = 'Drag Force [N]';
%save_figures(fig, './', ['mirror_predMap_'], fontsize, [7 6]);
