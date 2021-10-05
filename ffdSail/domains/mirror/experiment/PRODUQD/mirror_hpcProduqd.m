function mirror_hpcProduqd(encoding, nCases, startCase)
% Runs SAIL mirror optimization on cluster - [RUN THROUGH QSUB]
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: alexander.hagg@h-brs.de
% Oct 2017; Last revision: 23-Mar-2018

%------------- BEGIN CODE --------------
if nargin < 1
    encoding ='ffd'; nCases=5; startCase=1;
end

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
domainname = 'MIRROR'; % domainname = 'FOILFFD';
systemInit;
hpc = true;
if hpc; runFolder = '/home/ahagg2s/produqd'; else; runFolder = '/home/alex/produqd'; end
d = mirror_Domain('hpc', hpc, 'lowres', false, 'mirrorOnly', true, 'nCases', nCases, 'runFolder', runFolder);

% Algorithm hyperparameters
p = produqd;                                                % load default hyperparameters

% ADJUSTMENT FOR TESTING
p.qd.nGens              = 512;
p.qd.nInitialSamples    = 100;
p.qd.nAdditionalSamples = 10;
p.qd.nTotalSamples      = 400;
p.numIterations         = 1;

% Edit hyperparameters
p.qd.data.mapEval           = true;                         % produce intermediate prediction maps
p.qd.data.mapEvalMod        = p.qd.nAdditionalSamples;      % every n samples

% Set acquisition and prediction model parameters
for target=1:size(d.modelParsAcq,1)
    d.paramsAcq{target}     = feval(['params' d.modelParsAcq{target,2}], d.modelParsAcq{target,3:end});
    d.paramsPred{target}    = feval(['params' d.modelParsPre{target,2}], d.modelParsPre{target,3:end});
end
p.qd.data.predMapRes    = d.featureRes;

%% 1. INIT     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileINIT = ['produqd' domainname '_1_INIT.mat'];
if ~exist(fileINIT,'file')
    [observation, value] = initialSampling(d,p.qd.nInitialSamples);
    d.initialSampleSource = ['testSamples_100PE.mat'];
    save(d.initialSampleSource, 'observation', 'value');
    d.loadInitialSamples = true;
    save(fileINIT,'observation','value','p','d');
else
    load(fileINIT);
end

%% 2. SAIL     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileSAIL = ['produqd' domainname '_2_SAIL.mat'];
if ~exist(fileSAIL,'file')    
    runTime = tic;
    disp('Running PRODUQD first time to get samples and first SAIL run');
    p.numIterations         = 1;
    output = produqd(p,d);
    p.sailInput = output{end}; % Save for later PRODUQD runs
    disp(['Runtime: ' seconds2human(toc(runTime))]);
    save(fileSAIL,'output','p','d');
else
    load(fileSAIL);
    p.sailInput = output{end};
end

%% 3a. PRODUQD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.selectCriterion.type = 'prototypedist';
p.qd.nInitialSamples    = 400;
p.qd.nTotalSamples      = 700;
p.numIterations         = 2;

filePROD1 = ['produqd' domainname '_3a_PRODUQD.mat'];
if ~exist(filePROD1,'file')    
	disp('Running PRODUQD 1/2');
	p.selectCriterion.value = [1];
	runTime = tic;
	produqd_output{1} = produqd(p,d);
	disp(['Runtime: ' seconds2human(toc(runTime))]);
	save(filePROD1,'output','produqd_output','p','d');
else
    load(filePROD1);
end	

%% 3b. PRODUQD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filePROD2 = ['produqd' domainname '_3b_PRODUQD.mat'];
if ~exist(filePROD2,'file')    
	disp('Running PRODUQD 2/2');
	p.selectCriterion.value = [2];
	runTime = tic;
	produqd_output{2} = produqd(p,d);
	disp(['Runtime: ' seconds2human(toc(runTime))]);
	save(filePROD2,'output','produqd_output','p','d');
else
    load(filePROD2);
end	

%% 4. GT       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Getting ground truth of PRODUQD prediction maps');
d.nCases = nCases; % Make sure this is still correct... harhar
for run=1:length(produqd_output)
    %for it=1:length(produqd_output{run})
		it = 3;
        predMap = getGroundTruthMapArrayMirror(produqd_output{run}{it}.sail.predMap(end), d, 'dummy', false);
        produqd_output{run}{it}.sail.predMap(end).fitness_true = predMap.fitness_true;
        produqd_output{run}{it}.sail.predMap(end).dragForce_true = predMap.dragForce_true;
		save(['produqd' domainname '_4_GT.mat'],'output','produqd_output','p','d');
    %end
end

save(['produqd' domainname '_4_GT.mat'],'output','produqd_output','p','d');

disp('======== = DONE = ========')
end

%------------- END OF CODE --------------
