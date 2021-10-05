function FFDFOIL_hpcProduqd(encoding, nCases, startCase)
% Runs SAIL mirror optimization on cluster - [RUN THROUGH QSUB]
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: adam.gaier@h-brs.de
% Oct 2017; Last revision: 05-Oct-2017

%------------- BEGIN CODE --------------
if nargin < 1
    encoding ='ffd'; nCases=10;startCase=1;
end
%% Parallelization Settings
%parpool(12);

% Ensure Randomness of Randomness
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

% Create Temp Directory for Multithreading
tmpdir = getenv('TMPDIR');
if isempty(tmpdir);tmpdir='/tmp';end
myCluster.JobStorageLocation  = tmpdir;
myCluster.HasSharedFilesystem = true;

%% Add all files to path
addpath(genpath('~/produqd/'));
cd ~/produqd;

% Domain
domainname = 'FOILFFD';
systemInit;
hpc = true;
d = ffd_Domain('/scratch/ahagg2s/sailXFOIL');

% Algorithm hyperparameters
p = produqd;                                                % load default hyperparameters

% ADJUSTMENT FOR TESTING
p.qd.nGens              = 512;
p.qd.nInitialSamples    = 20;
p.qd.nAdditionalSamples = 10;
p.qd.nTotalSamples      = 120;
p.numIterations         = 3;

% Edit hyperparameters
p.qd.data.mapEval           = true;                         % produce intermediate prediction maps
p.qd.data.mapEvalMod        = 10;      % every n samples

% Set acquisition and prediction model parameters
for target=1:size(d.modelParsAcq,1)
    d.paramsAcq{target}     = feval(['params' d.modelParsAcq{target,2}], d.modelParsAcq{target,3:end});
    d.paramsPred{target}    = feval(['params' d.modelParsPre{target,2}], d.modelParsPre{target,3:end});
end
p.qd.data.predMapRes    = d.featureRes;

%% Run PRODUQD

runTime = tic;
disp('Running PRODUQD');
output = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
save(['produqd' domainname '.mat'],'output','p','d');

disp('Getting ground truth of PRODUQD prediction maps');
for it=1:length(output)
    predMap = getGroundTruthMapArrayMirror(output{it}.sail.predMap(end), d, 'dummy', false);
    output{it}.sail.predMap(end).fitness_true = predMap.fitness_true;
    output{it}.sail.predMap(end).dragForce_true = predMap.dragForce_true;
end

save(['produqd' domainname '_' int2str(randi(100000))  '.mat'],'output','p','d');

disp('======== = DONE = ========')
end

%------------- END OF CODE --------------
