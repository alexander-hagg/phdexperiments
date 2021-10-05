% Explore prototype-feature maps

%% Parallelization Settings
encoding ='ffd'; nCases=4;startCase=1;
%parpool(12);

% Ensure Randomness of Randomness
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

% Create Temp Directory for Multithreading
tmpdir = getenv('TMPDIR');
if isempty(tmpdir);tmpdir='/tmp';end
myCluster.JobStorageLocation  = tmpdir;
myCluster.HasSharedFilesystem = true;

% Add all files to path
addpath(genpath('~/produqd/'));
cd ~/produqd;

% Domain
domainname                      = 'FOILFFD';
systemInit;
hpc                             = false;
d                               = ffd_Domain('/scratch/ahagg2s/sailXFOIL');

% Algorithm hyperparameters
p                               = prodigi; % load default hyperparameters

% ADJUSTMENT FOR TESTING
p.qd.nGens                      = 512;
p.qd.nInitialSamples            = 100;
p.qd.nAdditionalSamples         = 10;
p.qd.nTotalSamples              = 1000;
p.numIterations                 = 1;
p.qd.constraints.threshold      = 0; %nan: no constraints, 0: hard constraints, > 0: soft constraints (* sigma)
% Edit hyperparameters
p.qd.data.mapEval               = true;    % produce intermediate prediction maps
p.qd.data.mapEvalMod            = 10;      % every n samples

% Set acquisition and prediction model parameters
for target=1:size(d.modelParsAcq,1)
    d.paramsAcq{target}         = feval(['params' d.modelParsAcq{target,2}], d.modelParsAcq{target,3:end});
    d.paramsPred{target}        = feval(['params' d.modelParsPre{target,2}], d.modelParsPre{target,3:end});
end
p.qd.data.predMapRes            = d.featureRes;

p.selectCriterion.value         = [1,2,3,4];

%% Run experiments
% Run PRODUQD with hard constraints
p.numIterations                 = 1;
runTime = tic;
disp('PRODUQD prior');
output = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
predMap{1} = output{1}.sail.predMap(end);
save(['produqd' domainname '.mat'],'output','p','d');


%% Configure prototype distance BC
d2 = d;
%output{1}.p.qd.constraintModel = trainConstraints(output{1}.optima,output{1}.latent);
%output{1}.p.qd.concept.allLabels = output{1}.conceptLabels
%output{1}.p.qd.oldconceptID = output{1}.conceptSelection.id;
d2.features.constraintModel = output{end}.p.qd.constraintModel;
selectID = output{end}.p.qd.concept.id;
labels = output{end}.p.qd.concept.allLabels;
d2.features.prototypes = output{end}.prototypesLatent(selectID,:);

d2.categorize                    = 'ffd_CategorizePrototypeDistance';
d2.featureMin                    = [0 0];
d2.featureMax                    = [1 1];
d2.featureLabels                 = {'d_{1,2}','d_{3,4}'};

% Create prediction map
disp(['Creating prediction map']);
qdcfg = output{1}.p.qd;
qdcfg.nGens = 1024;
qdcfg.constraints.threshold = nan;
qdcfg.data.predMapRes = [50 50];
[predMap{2}, percImproved(1,:)] = createPredictionMap(output{1}.sail.modelPred,qdcfg,d2,'featureRes',p.qd.data.predMapRes);
disp(['Percent improved in last iteration: ' num2str(percImproved(1,end))]);

qdcfg.constraints.threshold = 0;
[predMap{3}, percImproved(2,:)] = createPredictionMap(output{1}.sail.modelPred,qdcfg,d2,'featureRes',p.qd.data.predMapRes);
disp(['Percent improved in last iteration: ' num2str(percImproved(2,end))]);

%qdcfg.constraints.threshold = 2;
%[predMap{4}, percImproved(3,:)] = createPredictionMap(output{1}.sail.modelPred,qdcfg,d2,'featureRes',p.qd.data.predMapRes);
%disp(['Percent improved in last iteration: ' num2str(percImproved(3,end))]);

disp('======== = DONE = ========')

%% Get Ground Truth
disp('===============        GT       ================');
for i=1:length(predMap)
    if i==1
        dom = d;
    else
        dom = d2;
    end
    tMap = getGroundTruthMapArray(predMap{i}, dom, 'dummy', false);
    predMap{i}.fitness_true = tMap.fitness_true;
    predMap{i}.cD_true = tMap.cD_true;
    predMap{i}.cL_true = tMap.cL_true;
end
%% Analysis


fig(1) = figure(1);
viewMap(predMap{1}.fitness_true,d)
caxis([-5.5 -2]);title('Base Feature Space');
fig(2) = figure(2);
viewMap(predMap{2}.fitness_true,d2)
caxis([-5.5 -2]);title('Prototype Distance Space, No Selection');
fig(3) = figure(3);
viewMap(predMap{3}.fitness_true,d2)
caxis([-5.5 -2]);title('Prototype Distance Space, Selection');
% fig(4) = figure(4);
% viewMap(predMap{4}.fitness_true,d2)
% caxis([-5.5 -2]);title('Prototype Distance Space, Soft Constr.');

save_figures(fig, '.', 'PrototypeDistanceMaps', 14, [7 7]);

%%

figure(5);
errors1 = abs(predMap{1}.fitness_true-predMap{1}.fitness);
errors2 = abs(predMap{2}.fitness_true-predMap{2}.fitness);
errors3 = abs(predMap{3}.fitness_true-predMap{3}.fitness);
errors4 = abs(predMap{4}.fitness_true-predMap{4}.fitness);
boxplot([errors1(:) errors2(:) errors3(:) errors4(:)]);
title('Abs Error');
grid on;
axis([0.5 4.5 -0.1 0.5]);

%%
% figure(3);
% subplot(2,3,1);hold off;
% foil = output2{1}.prototypes(output2{1}.conceptSelection.id(1),:);
% visFoil(foil');
% title('Prototype 1');
% subplot(2,3,2);hold off;
% foil = output2{1}.prototypes(output2{1}.conceptSelection.id(2),:);
% visFoil(foil');
% title('Prototype 2');
% subplot(2,3,3);hold off;
% foil = output2{1}.prototypes(output2{1}.conceptSelection.id(3),:);
% visFoil(foil');
% title('Prototype 3');
%
% subplot(2,3,4);hold off;
% visFoil(squeeze(output2{2}.sail.predMap(end).genes(1,1,:)));
% features = ffd_CategorizePrototypeDistance(squeeze(output2{2}.sail.predMap(end).genes(1,1,:))',d2);
% title([num2str(features)]);
% subplot(2,3,5);hold off;
% visFoil(squeeze(output2{2}.sail.predMap(end).genes(20,25,:)));
% features = ffd_CategorizePrototypeDistance(squeeze(output2{2}.sail.predMap(end).genes(20,25,:))',d2);
% title([num2str(features)]);
% subplot(2,3,6);hold off;
% visFoil(squeeze(output2{2}.sail.predMap(end).genes(25,21,:)));
% features = ffd_CategorizePrototypeDistance(squeeze(output2{2}.sail.predMap(end).genes(25,21,:))',d2);
% title([num2str(features)]);

%%

