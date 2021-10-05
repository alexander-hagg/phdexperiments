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
p.qd.nGens              = 512;
p.qd.nInitialSamples    = 100;
p.qd.nAdditionalSamples = 10;
p.qd.nTotalSamples      = 300;
p.numIterations         = 3;
p.qd.constraints.threshold   = nan; %nan: no constraints, 0: hard constraints, > 0: soft constraints (* sigma)
% Edit hyperparameters
p.qd.data.mapEval       = true;    % produce intermediate prediction maps
p.qd.data.mapEvalMod    = 10;      % every n samples

% Set acquisition and prediction model parameters
for target=1:size(d.modelParsAcq,1)
    d.paramsAcq{target}     = feval(['params' d.modelParsAcq{target,2}], d.modelParsAcq{target,3:end});
    d.paramsPred{target}    = feval(['params' d.modelParsPre{target,2}], d.modelParsPre{target,3:end});
end
p.qd.data.predMapRes    = d.featureRes;

p.selectCriterion.value = 5;
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



%% Run good old fashioned PRODUQD
runTime = tic;
disp('Running PRODUQD without constraints');
p.qd.constraints.threshold = nan;
output{1} = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
save(['produqd' domainname '.mat'],'output','p','d');

% Run PRODUQD with hard constraints
runTime = tic;
disp('Running PRODUQD with hard constraints');
p.qd.constraints.threshold = 0;
output{2} = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
save(['produqd' domainname '.mat'],'output','p','d');

% Run PRODUQD with soft constraints
runTime = tic;
disp('Running PRODUQD with soft constraints');
p.qd.constraints.threshold = 0.5;
output{3} = produqd(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
save(['produqd' domainname '.mat'],'output','p','d');

% Run PRODUQD with soft constraints
runTime = tic;
disp('Running PRODUQD with soft constraints');
p.qd.constraints.threshold = 1;
output{4} = prodigi(p,d);
disp(['Runtime: ' seconds2human(toc(runTime))]);
save(['produqd' domainname '.mat'],'output','p','d');

% % Run PRODUQD with soft constraints
% runTime = tic;
% disp('Running PRODUQD with soft constraints');
% p.qd.constraints.threshold = 1.5;
% output{5} = produqd(p,d);
% disp(['Runtime: ' seconds2human(toc(runTime))]);
% save(['produqd' domainname '.mat'],'output','p','d');


%%
% disp('Getting ground truth of PRODUQD prediction maps');
% for it=1:length(output)
%     predMap = getGroundTruthMapArrayMirror(output{it}.sail.predMap(end), d, 'dummy', false);
%     output{it}.sail.predMap(end).fitness_true = predMap.fitness_true;
%     output{it}.sail.predMap(end).dragForce_true = predMap.dragForce_true;
% end
%
% save(['produqd' domainname '_' int2str(randi(100000))  '.mat'],'output','p','d');

disp('======== = DONE = ========')

%% Analysis
dists = []; runIDs = [];
for i=1:length(output)
    id = output{i}{1}.conceptSelection.id;
    prototype = output{i}{1}.prototypes(id,:);
    
    for j=1:length(output{i})
        if i>1&&j==1; continue; end
        newdists = pdist2(output{i}{j}.optima,prototype);
        dists = [dists; newdists];
        if j==1
            runIDs = [runIDs;((i-1)*2+1)*ones(length(newdists),1)];
        else
            runIDs = [runIDs;((i-1)*2+2)*ones(length(newdists),1)];
        end
    end
end

fig(1) = figure(1);
boxplot(dists,runIDs,'Notch','on');
ax = gca;
ylabel('\delta Parameter Space');
grid on;
ax.XTickLabel = {'N0', 'N1+', 'H1+','S0.5|1+','S1.0|1+'};

newdists = []; dists = []; runIDs = [];
for i=1:length(output)
    id = output{i}{1}.conceptSelection.id
    prototype = output{i}{1}.prototypesLatent(id,:);
    for j=1:length(output{i})
        if i>1&&j==1; continue; end
        newdists = pdist2(output{i}{j}.latent,prototype);
        dists = [dists; newdists];
        if j==1
            runIDs = [runIDs;((i-1)*2+1)*ones(length(newdists),1)];
        else
            runIDs = [runIDs;((i-1)*2+2)*ones(length(newdists),1)];
        end
    end
end

fig(2) = figure(2);
boxplot(dists,runIDs,'Notch','on');
ax = gca;
ylabel('\delta Similarity Space');
ax.XTickLabel = {'N0', 'N1+', 'H1+','S0.5|1+','S1.0|1+'};
grid on;

dists = []; runIDs = [];
for i=1:length(output)
    id = output{i}{1}.conceptSelection.id;
    prototype = output{i}{1}.prototypes(id,:);
    
    for j=1:length(output{i})
        if i>1&&j==1; continue; end
        newdists = pdist2(output{i}{j}.optima,prototype);
        dists = [dists; newdists];
        runIDs = [runIDs;((i-1)*length(output{i})+j)*ones(length(newdists),1)];
    end
end

fig(3) = figure(3);
boxplot(dists,runIDs,'Notch','on');
ax = gca;
ylabel('\delta Parameter Space');
grid on;
ax.XTickLabel = {'N0', 'N1', 'N2', 'H1', 'H2', 'S0.5|1', 'S0.5|2', 'S1.0|1', 'S1.0|2'};

newdists = []; dists = []; runIDs = [];
for i=1:length(output)
    id = output{i}{1}.conceptSelection.id
    prototype = output{i}{1}.prototypesLatent(id,:);
    for j=1:length(output{i})
        if i>1&&j==1; continue; end
        newdists = pdist2(output{i}{j}.latent,prototype);
        dists = [dists; newdists];
        runIDs = [runIDs;((i-1)*length(output{i})+j)*ones(length(newdists),1)];
    end
end

fig(4) = figure(4);
boxplot(dists,runIDs,'Notch','on');
ax = gca;
ylabel('\delta Similarity Space');
ax.XTickLabel = {'N0', 'N1', 'N2', 'H1', 'H2', 'S0.5|1', 'S0.5|2', 'S1.0|1', 'S1.0|2'};
grid on;

save_figures(fig, '.', 'constraints', 14, [10 4]);

%%
cmap = hsv(3);
clear fig;
for i=1:length(output)
    fig(i) = figure(i);hold off;
    for j=1:length(output{i})
        l = output{i}{j}.latent;
        h(j) = scatter(l(:,1),l(:,2),128, cmap(j,:), 'filled'); hold on;
        h(j).MarkerFaceAlpha = 0.04;
    end
    
    for j=1:length(output{i})
        lp = output{i}{j}.prototypesLatent;
        h(j) = scatter(lp(:,1),lp(:,2),32, cmap(j,:), 'filled'); hold on;
        h(j).MarkerEdgeColor = 'w';
    end
    
    for j=1:length(output{i})
        lp = output{i}{j}.prototypesLatent;
        sel = output{i}{j}.conceptSelection.id;
        d = scatter(lp(sel,1),lp(sel,2),256, cmap(j,:), 'd', 'filled'); hold on;
        d.MarkerEdgeColor = 'w';
        t = text(lp(sel,1)+1,lp(sel,2)+1,['P' int2str(j)]);
        %t.Color = 'k';
    end
    legend([h(1) h(2) h(3)], '1', '2''', '3''');
    axis equal;
    axis([-40 40 -20 20]);
end

%%

clear fig;close all
for i=1:length(output)
    fig(i) = figure(i);hold off;
    for j=1:length(output{i})
        sel = output{i}{j}.conceptSelection.id;
        oo = output{i}{j}.prototypes(sel,:);
        %base = loadBaseAirfoil(d.express, 10, '.');
        %plot(base.foil(1,:),base.foil(2,:), 'Color', [0 0 0], 'LineWidth',2);hold on;
        %hold on;
        %drawnow;
        for g=1:size(oo,1)
            foil = ffdRaeY(oo(g,:));
            plot1(j) = plot(foil(1,:),foil(2,:), 'LineWidth',4);hold on;
            %plot1.Color(4) = 0.5;
        end
        axis equal;
        axis tight;
        axis([0 1 -0.15 0.15]);
        l = legend('1', '2', '3');
        l.FontSize = 14;
    end
end
save_figures(fig, '.', 'constraints_foilprototypes', 14, [10 4]);

end

%------------- END OF CODE --------------
