function FFDFOIL_inclass_acquisition(encoding, nCases, startCase)
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
    encoding ='ffd'; nCases=4;startCase=1;
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

% Add all files to path
addpath(genpath('~/produqd/'));
cd ~/produqd;

% Domain
domainname = 'FOILFFD';
systemInit;
hpc = false;
d = ffd_Domain('/scratch/ahagg2s/sailXFOIL');

% Algorithm hyperparameters
p = produqd;                                                % load default hyperparameters

% ADJUSTMENT FOR TESTING
p.qd.nGens                  = 512;
p.qd.nInitialSamples        = 100;
p.qd.nAdditionalSamples     = 10;
p.qd.nTotalSamples          = 300;
p.numIterations             = 2;
p.qd.constraints.threshold  = nan; %nan: no constraints, 0: hard constraints, > 0: soft constraints (* sigma)
p.qd.data.mapEval           = true;    % produce intermediate prediction maps
p.qd.data.mapEvalMod        = 10;      % every n samples

% Set acquisition and prediction model parameters
for target=1:size(d.modelParsAcq,1)
    d.paramsAcq{target}     = feval(['params' d.modelParsAcq{target,2}], d.modelParsAcq{target,3:end});
    d.paramsPred{target}    = feval(['params' d.modelParsPre{target,2}], d.modelParsPre{target,3:end});
end
p.qd.data.predMapRes    = d.featureRes;

p.qd.display.illu = false;
p.qd.display.illuMod = 25;

p.selectCriterion.value = 3;
%% Run experiments
% Get first SAIL run for comparabiltiy
disp('===============   SAIL PRIOR    ================');
runTime = tic;
p.sailInput = sail(p.qd,d);
disp(['Runtime SAIL prior: ' seconds2human(toc(runTime))]);

% Determine similarity space
p.simspaceInput.optima = reshape(p.sailInput.predMap(end).genes,prod(p.qd.data.predMapRes),d.dof);
p.simspaceInput.optima(any(isnan(p.simspaceInput.optima')),:) = [];
% Extract concepts
[p.simspaceInput.conceptLabels, p.simspaceInput.latent, ~, ~, ~, p.simspaceInput.concepts] = ...
    dimReducedClustering( p.simspaceInput.optima, 'tSNE', 2);

% Extract prototypes, excluding non-assigned concept ID '0'
for ii=1:max(p.simspaceInput.concepts.uniqid)
    samples = p.simspaceInput.optima(p.simspaceInput.conceptLabels==ii,:);
    latent = p.simspaceInput.latent(p.simspaceInput.conceptLabels==ii,:);
    [~,p.simspaceInput.prototypes(ii,:),~,~,protIDs] = kmedoids(samples,1);
    p.simspaceInput.prototypesLatent(ii,:) = latent(protIDs,:);
end

% Run good old fashioned PRODUQD
runTime = tic;
disp('Running PRODUQD without constraints');
p.qd.constraints.threshold = nan;
output{1} = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
save(['produqd' domainname '.mat'],'output','p','d');

% Run PRODUQD with hard constraints w/o acquisition outclass penalty
runTime = tic;
disp('Running PRODUQD with hard constraints, no outclass penalty');
d.fitnessPenalty_OutClass = false;
p.qd.constraints.threshold = 0;
output{2} = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
save(['produqd' domainname '.mat'],'output','p','d');

% Run PRODUQD with hard constraints with acquisition outclass penalty
runTime = tic;
disp('Running PRODUQD with hard constraints, incl. outclass penalty');
d.fitnessPenalty_OutClass = true;
p.qd.constraints.threshold = 0;
output{3} = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
save(['produqd' domainname '.mat'],'output','p','d');

%
disp('Getting ground truth of PRODUQD prediction maps');
for o=1:length(output)
    for it=1:length(output{o})
        predMap = getGroundTruthMapArray(output{o}{it}.sail.predMap(end), d, 'dummy', false);
        output{o}{it}.sail.predMap(end).cD_true = predMap.cD_true;
        output{o}{it}.sail.predMap(end).cL_true = predMap.cL_true;
        output{o}{it}.sail.predMap(end).fitness_true = predMap.fitness_true;        
    end
end

save(['produqd' domainname '_FINISHED.mat'],'output','p','d');

disp('======== = DONE = ========')

%% Visualization
close all; clear fig; set(0,'DefaultFigureWindowStyle','docked');

for o=1:length(output)
    for it=1:length(output{o})
        fig((o-1)*4+it) = figure((o-1)*4+it);hold off;
        viewMap(output{o}{it}.sail.predMap(end).fitness_true,d);
        caxis([-5 -3]);
    end
end



%%
distances = [];distancesNOT = [];grps = [];
parameterVariation=[];grpsParVar = [];
grpID = 0;
clear pts;
for it=1:length(output)
    disp('..')
    if it==1
        model = trainConstraints(output{1}{1}.optima,output{1}{1}.latent);
        start = 1;
    else
        model = output{it}{1}.p.qd.constraints.model;
        start = 2;
    end
    for iter=start:length(output{it})
        protoID = output{it}{1}.conceptSelection.id;
        centers = output{it}{1}.prototypesLatent; center = centers(protoID,:);
        %XXX = reshape(output{it}{iter}.sail.predMap(end).genes,625,10);
        if iter==1
            XXX = output{it}{iter}.sail.model{1}.trainInput;
        else
            XXX = output{it}{iter}.sail.model{1}.trainInput(end-200:end,:);
        end
        [m_1, s2_1] = gp(model.x.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
            model.trainInput, model.trainOutput(:,1), XXX);
        [m_2, s2_2] = gp(model.y.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
            model.trainInput, model.trainOutput(:,2), XXX);
        pts = [m_1,m_2];
    
        distances = [distances pdist2(center,pts)];
        distancesNOT = [distancesNOT pdist2(center,pts)];
        parameterVariation = [parameterVariation std(XXX)];
        grps = [grps grpID*ones(1,length(pdist2(center,pts)))];
        grpsParVar = [grpsParVar grpID*ones(1,length(std(XXX)))];
        grpID = grpID + 1;
    end
end


close all; clear fig; set(0,'DefaultFigureWindowStyle','default');

clear fig;
fig(1) = figure(1);
h = boxplot(distances,grps,'PlotStyle','compact' , 'Colors', 'k');grid on;
ylabel('Distance to Prototype');
set(h(5,:),'Visible','off')
ax = gca;
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'1st', 'No', 'Hard', 'InclassAcq'};
title('Samples');
drawnow;

save_figures(fig, './', ['inclass_acquisition_boxplots'], 18, [6 4]);
%% Boxplot fitness and distance prediction map


fitValues = [];grps = [];distances = [];
grpID = 0;
clear pts;
for it=1:length(output)
    disp('..')
    if it==1
        start = 1;
    else
        start = 2;
    end
    for iter=start:length(output{it})
        % Fitness values
        fits = output{it}{iter}.sail.predMap(end).fitness_true(:);
        fitValues = [fitValues; fits(~isnan(fits))];
        
        % Distances
        XXX = reshape(output{it}{iter}.sail.predMap(end).genes,625,10);
        XXX(isnan(fits),:) = [];
        XXX(any(isnan(XXX')),:) = [];
        [m_1, s2_1] = gp(model.x.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
            model.trainInput, model.trainOutput(:,1), XXX);
        [m_2, s2_2] = gp(model.y.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
            model.trainInput, model.trainOutput(:,2), XXX);
        pts = [m_1,m_2];
    
        distances = [distances pdist2(center,pts)];
        
        % Grouping
        grps = [grps grpID*ones(1,length(fits(~isnan(fits))))];
        grpID = grpID + 1;
    end
end

close all; clear fig; set(0,'DefaultFigureWindowStyle','default');
fig(1) = figure(1);
h = boxplot(fitValues,grps,'PlotStyle','compact' , 'Colors', 'k');grid on;
ylabel('True Fitness Values');
set(h(5,:),'Visible','off')
ax = gca;
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'1st', 'No', 'Hard', 'InclassFit'};
axis([0.5 4.5 -5 -3.5]);
title('Prediction Map');
drawnow;


fig(2) = figure(2);
h = boxplot(distances,grps,'PlotStyle','compact' , 'Colors', 'k');grid on;
ylabel('Distance to Prototype');
set(h(5,:),'Visible','off')
ax = gca;
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'1st', 'No', 'Hard', 'InclassFit'};
axis([0.5 4.5 0 45]);
title('Prediction Map');

save_figures(fig, './', ['inclass_predictionMap_boxplots'], 18, [6 4]);


%% Latent space

for it=1:length(output)
    if it==1
        model = trainConstraints(output{1}{1}.optima,output{1}{1}.latent);
    else
        model = output{it}{1}.p.qd.constraints.model;
    end
    smp = output{it}{end}.sail.model{1}.trainInput;
    opt = reshape(output{it}{iter}.sail.predMap(end).genes,625,10);
    opt(any(isnan(opt')),:) = [];
    samples{it} = [smp; opt];
    grpsLatent{it} = [zeros(size(smp,1),1); ones(size(opt,1),1)];
    
    [m_1, s2_1] = gp(model.x.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
            model.trainInput, model.trainOutput(:,1), samples{it});
    [m_2, s2_2] = gp(model.y.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
            model.trainInput, model.trainOutput(:,2), samples{it});
    newLatent{it} = [m_1,m_2];
end

titles = {'Seeding Only', 'Hard Constraints', 'Hard Constraints with constrained acquisition'};
for i=1:length(output)
    fig(i) = figure(i);hold off;
    cmap = parula(2);
    
    pts = newLatent{i}(grpsLatent{i}==0,:);
    s1a = scatter(pts(1:300,1),pts(1:300,2),16,[0.7 0.7 0.7],'filled');hold on;
    s1b = scatter(pts(301:end,1),pts(301:end,2),16,[0 0 0 ],'filled');hold on;
    pts = newLatent{i}(grpsLatent{i}==1,:);
    s2 = scatter(pts(:,1),pts(:,2),16,[1 0.3 0.3],'filled');hold on;
    title(titles{i});
    legend([s1a s1b s2],'Samples Iter 1', 'Samples Iter 2', 'Optima Iter 2');
    ax = gca;
    ax.Visible = 'off';
end

save_figures(fig, './', ['inclass_acquisition'], 12, [6 4]);


%------------- END OF CODE --------------
end
