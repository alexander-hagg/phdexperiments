%PRODUQD vs SAIL comparison -
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
domainname = 'FOILFFD';
systemInit;
d = ffd_Domain('/scratch/ahagg2s/sailXFOIL');

numPRODUQDRuns = 10;

% Algorithm hyperparameters
p = produqd;                                                % load default hyperparameters

% ADJUSTMENT FOR TESTING
p.qd.nGens              = 512;
p.qd.nInitialSamples    = 50;
p.qd.nAdditionalSamples = 10;
p.qd.nTotalSamples      = 200;
p.numIterations         = 3;

% Edit hyperparameters
p.qd.data.mapEval           = true;                         % produce intermediate prediction maps
p.qd.data.mapEvalMod        = p.qd.nAdditionalSamples;      % every n samples

% Set acquisition and prediction model parameters
for target=1:size(d.modelParsAcq,1)
    d.paramsAcq{target}     = feval(['params' d.modelParsAcq{target,2}], d.modelParsAcq{target,3:end});
    d.paramsPred{target}    = feval(['params' d.modelParsPre{target,2}], d.modelParsPre{target,3:end});
end
p.qd.data.predMapRes    = d.featureRes;

%% Get samples
[observation, value] = initialSampling(d,p.qd.nInitialSamples);
d.initialSampleSource = 'initSamples.mat';
save(d.initialSampleSource, 'observation', 'value');
d.loadInitialSamples = true;

%% Run SAIL first time
disp('===============   SAIL PRIOR    ================');
runTime = tic;
qd = p.qd;
p.sailInput = sail(qd,d);
disp(['Runtime SAIL prior: ' seconds2human(toc(runTime))]);
p.qd.nInitialSamples    = p.qd.nTotalSamples;

% Save samples
observation = [p.sailInput.tModelsAcq{end,1}.trainInput];
value = nan(size(observation,1),size(p.sailInput.tModelsAcq,2));
for i = 1:size(p.sailInput.tModelsAcq,2)
    value(:,i) = p.sailInput.tModelsAcq{end,i}.trainOutput;
end
d.initialSampleSource = 'priorSamples.mat';
save(d.initialSampleSource, 'observation', 'value');

%% Get similarity space, classes and prototypes

optima = reshape(p.sailInput.predMap(end).genes,prod(p.qd.data.predMapRes),d.dof);
optima(any(isnan(optima')),:) = [];

% Extract concepts
[conceptLabels, latent, ~, ~, ~, concepts] = ...
    dimReducedClustering( optima, 'tSNE', 2);

% Extract prototypes, excluding non-assigned concept ID '0'
for ii=1:max(concepts.uniqid)
    samplesLatent = latent(conceptLabels==ii,:);
    samples = optima(conceptLabels==ii,:);
    [~,prototypesLatent(ii,:),~,~,protIDs] = kmedoids(samplesLatent,1);
    prototypes(ii,:) = samples(protIDs,:);
end


p.simspaceInput.optima = optima;
p.simspaceInput.conceptLabels = conceptLabels;
p.simspaceInput.latent = latent;
p.simspaceInput.concepts = concepts;
p.simspaceInput.prototypes = prototypes;
p.simspaceInput.prototypesLatent = prototypesLatent;

%% Continue SAIL for comparison
disp('===============   SAIL CONT.    ================');
runTime = tic;
qd = p.qd; qd.nTotalSamples = p.qd.nTotalSamples * (p.numIterations); % SAIL gets as many PE as PRODUQD
disp(['Total samples: ' int2str(qd.nTotalSamples)]);
output{1}.data{1}.sail = sail(qd,d);
disp(['Runtime SAIL: ' seconds2human(toc(runTime))]);

%% Run PRODUQD
disp('===============     PRODUQD     ================');
runTime = tic;
p.qd.nTotalSamples = p.qd.nTotalSamples * 2;
disp(['Total samples: ' int2str(p.qd.nTotalSamples) ' and ' int2str(p.numIterations) ' iterations']);
for critVal=1:numPRODUQDRuns
    p_crit = p; d_crit = d;
    p_crit.selectCriterion.value = critVal;
    output{1+critVal}.data = produqd(p_crit,d_crit);
    disp(['Runtime PRODUQD run ' int2str(critVal) ': ' seconds2human(toc(runTime))]);
end
%% Get Ground Truth
disp('===============        GT       ================');
predMap = getGroundTruthMapArray(output{1}.data{end}.sail.predMap(end), d, 'dummy', false);
output{1}.data{end}.sail.predMap(end).fitness_true = predMap.fitness_true;
output{1}.data{end}.sail.predMap(end).cD_true = predMap.cD_true;
output{1}.data{end}.sail.predMap(end).cL_true = predMap.cL_true;
for it=2:length(output)
    for i=2:3
        predMap = getGroundTruthMapArray(output{it}.data{i}.sail.predMap(end), d, 'dummy', false);
        output{it}.data{i}.sail.predMap(end).fitness_true = predMap.fitness_true;
        output{it}.data{i}.sail.predMap(end).cD_true = predMap.cD_true;
        output{it}.data{i}.sail.predMap(end).cL_true = predMap.cL_true;
    end
end

%% Save results
save(['SAILvsPRODUQD_' int2str(randi(10000))],'output','p','d');

%% Visualization
close all; clear fig; set(0,'DefaultFigureWindowStyle','default');

for i=2:numPRODUQDRuns+1
    fig(i-1) = figure(i-1);hold off;
    viewMap(output{1}.data{end}.sail.predMap(end).fitness_true - output{i}.data{end}.sail.predMap(end).fitness_true,d);
    caxis([-0.3 0.3]);
    cmap = parula(2);
    %cmap(2,:) = cmap(1,:);
    %cmap(4,:) = cmap(5,:);
    colormap(cmap);
end

save_figures(fig, './', ['sail_vs_produqd_FIT_'], 40, [9 7]);

%% Maps (Distance)
%close all;
clear fig; set(0,'DefaultFigureWindowStyle','default');

dists{1} = [];
dists{2} = [];
for it=2:length(output)
    disp('..')
    centersLatent = output{it}.data{1}.prototypesLatent;
    centers = output{it}.data{1}.prototypes;
    protoID = output{it}.data{1}.conceptSelection.id;
    centerLatent = centersLatent(protoID,:);
    center = centers(protoID,:);
    
    %XXX = reshape(output{1}.data{end}.sail.predMap(end).genes,625,10);
    XXX = reshape(output{it}.data{1}.sail.predMap(end).genes,625,10);
    dists{1} = pdist2(center,XXX);
    XXX = reshape(output{it}.data{end}.sail.predMap(end).genes,625,10);
    dists{2} = pdist2(center,XXX);
    
    % Get prototype location in feature map
    [~, index]=ismember(center,reshape(output{it}.data{1}.sail.predMap(end).genes,625,10),'rows');
    [protOrgX,protOrgY] = ind2sub(size(output{1}.data{end}.sail.predMap(end).fitness),index);
    
    distMap = reshape(dists{2},25,25);
    distMap(protOrgX,protOrgY) = 0;
    
    fig(it-1) = figure(it-1);
    [figHandle, imageHandle, cHandle] = viewMap(distMap,d);
    cHandle.Label.String = 'Distance to prototype';caxis([0 1.5]);colormap(parula(5));
    
end
for i=1:length(fig)
    figure(i)
    name = ['sail_vs_produqd_DIST_'];
    %save2pdf([name int2str(i) '.pdf'],fig(i),600,[7 6]);
    %print([name int2str(i)],'-depsc');
    saveas(gcf,['ffd_distmap_' int2str(i) '.eps'],'epsc')
end

%% Boxplots (distance, fitness, prediction error)

dists{1} = [];
dists{2} = [];
dists{3} = [];
for it=2:length(output)
    disp('..')
    %centersLatent = output{it}.data{1}.prototypesLatent;
    %centerLatent = centersLatent(protoID,:);
    
    for iter=1:length(output{it}.data)
        centers = output{it}.data{iter}.prototypes;
        protoID = output{it}.data{iter}.conceptSelection.id;
        center = centers(protoID,:);
        
        XXX = reshape(output{it}.data{iter}.sail.predMap(end).genes,625,10);
        dists{iter} = [dists{iter} pdist2(center,XXX)];
    end
    
    %XXX = output{it}.data{1}.latent;
    %dists{1} = [dists{1}, pdist2(centerLatent,XXX)];
    %XXX = output{it}.data{end}.latent;
    %dists{2} = [dists{2} pdist2(centerLatent,XXX)];
end

clear fig;
fig(1) = figure(1);
distances = [dists{1}(:);dists{2}(:);dists{3}(:)];
grps = [zeros(1,length(dists{1})),ones(1,length(dists{2})),2*ones(1,length(dists{3}))];
h = boxplot(distances,grps,'PlotStyle','compact' , 'Colors', 'k');grid on;
%axis([0.5 2.5 -0.01 2.2]);
ax = gca;
ax.XTick = [1 2 3]; ax.XTickLabel = {'SAIL', 'P 1', 'P 2'};
ylabel('Distance to Prototype');
set(h(5,:),'Visible','off')
drawnow;
clear fits
%
fitsTrue{1} = output{1}.data{end}.sail.predMap(end).fitness_true(:);
fits{1} = output{1}.data{end}.sail.predMap(end).fitness(:);
fitsTrue{2} = [];fitsTrue{3} = [];
fits{2} = [];fits{3} = [];
for it=2:length(output)
    fitsTrue{2} = [fitsTrue{2}; output{it}.data{2}.sail.predMap(end).fitness_true(:)];
    fitsTrue{3} = [fitsTrue{3}; output{it}.data{3}.sail.predMap(end).fitness_true(:)];
    fits{2} = [fits{2}; output{it}.data{2}.sail.predMap(end).fitness(:)];
    fits{3} = [fits{3}; output{it}.data{3}.sail.predMap(end).fitness(:)];
end

%
fig(2) = figure(2);
C = [fitsTrue{1}(:);fitsTrue{2}(:);fitsTrue{3}(:)];
grp = [zeros(1,length(fitsTrue{1}(:))),ones(1,length(fitsTrue{2}(:))), 2*ones(1,length(fitsTrue{2}(:)))];
h = boxplot(C,grp,'PlotStyle','compact', 'Colors', 'k');grid on;
axis([0.5 3.5 -5 -3.6]);
ax = gca;
ax.XTick = [1 2 3]; ax.XTickLabel = {'SAIL', 'P 1', 'P 2'};
ylabel('\leftarrow Fitness Score');
set(h(5,:),'Visible','off')
drawnow;
%
fig(3) = figure(3);
err{1} = abs(fits{1}-fitsTrue{1});
err{2} = abs(fits{2}-fitsTrue{2});
err{3} = abs(fits{3}-fitsTrue{3});

err{1}(isnan(err{1})) = [];
err{2}(isnan(err{2})) = [];
err{3}(isnan(err{3})) = [];

grp = [zeros(1,length(err{1})),ones(1,length(err{2})),2*ones(1,length(err{3}))];
h = boxplot([err{1};err{2};err{3}],grp,'PlotStyle','compact','Colors', 'k');grid on;
axis([0.5 3.5 -0.01 0.3]);
ax = gca;
ax.XTick = [1 2 3]; ax.XTickLabel = {'SAIL', 'P 1', 'P 2'};
ylabel('Prediction Error');
set(h(5,:),'Visible','off')
drawnow;

%% Boxplots inner variance PRODUQD vs variance SAIL
genes = reshape(output{1}.data{1}.sail.predMap(end).genes,25*25,10);
stdGenes{1} = nanstd(genes);
stdGenes{2} = [];stdGenes{3} = [];
for i=2:numPRODUQDRuns
    for iter=2:3
        genes = reshape(output{i}.data{iter}.sail.predMap(end).genes,25*25,10);
        labels = output{i}.data{iter}.conceptLabels;
        uniqLabels = unique(labels);
        uniqLabels(1) = [];
        for j=1:max(uniqLabels)
            stdGenes{iter}(end+1,:) = nanstd(genes(labels==j,:));
        end
    end
end

fig(4) = figure(4);
C = [stdGenes{1}(:); stdGenes{2}(:); stdGenes{3}(:)];
grp = [zeros(1,length(stdGenes{1}(:))) ones(1,length(stdGenes{2}(:))) 2*ones(1,length(stdGenes{3}(:)))];
h = boxplot(C,grp,'PlotStyle','compact', 'Colors', 'k');grid on;
axis([0.5 3.5 -0.01 0.4]);
ax = gca;
ax.XTick = [1 2 3]; ax.XTickLabel = {'SAIL', 'P 1', 'P 2'};
ylabel('\sigma Design Parameters');
set(h(5,:),'Visible','off')
drawnow;

save_figures(fig, './', ['sail_vs_produqd_boxplots_'], 20, [4 5]);

%% Latent space
for i=2:numPRODUQDRuns
    optima{i} = [];
    for j=1:p.numIterations
        optima{i} = [optima{i}; output{i}.data{j}.optima];
    end
    [~, newLatent{i}] = dimReducedClustering( optima{i}, 'tSNE', 2);
end
%%
clear grp
for i=2:numPRODUQDRuns
    grp{i} = [];
    for j=1:p.numIterations
        grp{i} = [grp{i}, (j-1)*ones(1,size(output{i}.data{j}.optima,1))];
    end
end
%%
for i=2:numPRODUQDRuns
    fig(i-1) = figure(i-1);hold off;
    cmap = parula(2);
    
    
    
    s2 = scatter(newLatent{i}(grp{i}==1,1),newLatent{i}(grp{i}==1,2),32,cmap(1,:),'d','filled');hold on;
    s3 = scatter(newLatent{i}(grp{i}==2,1),newLatent{i}(grp{i}==2,2),32,cmap(2,:),'^','filled');hold on;
    s1 = scatter(newLatent{i}(grp{i}==0,1),newLatent{i}(grp{i}==0,2),16,[0.7 0.7 0.7],'filled');hold on;
    
    labels1 = output{i}.data{1}.conceptLabels;
    labels2 = output{i}.data{2}.conceptLabels;
    labels3 = output{i}.data{3}.conceptLabels;
    sel1 = [[labels1==output{i}.data{1}.conceptSelection.id]; logical(zeros(length(labels2)+length(labels3),1))];
    s4 = scatter(newLatent{i}(logical(sel1),1),newLatent{i}(logical(sel1),2),32,[1 1 1],'r','filled');
    
    
    %     clear selIDs;
    %     labels1 = output{i}.data{1}.conceptLabels;
    %     labels2 = output{i}.data{2}.conceptLabels;
    %     labels3 = output{i}.data{3}.conceptLabels;
    %     sel1 = [[labels1==output{i}.data{1}.conceptSelection.id]; logical(zeros(length(labels2)+length(labels3),1))];
    %     sel2 = [logical(zeros(length(labels1),1)); [labels2==output{i}.data{2}.conceptSelection.id]; logical(zeros(length(labels3),1))];
    %     sel3 = [logical(zeros(length(labels1)+length(labels2),1)); [labels3==output{i}.data{3}.conceptSelection.id]];
    %
    %     cmap = parula(3);
    %     s1 = scatter(newLatent{i}(logical(sel1),1),newLatent{i}(logical(sel1),2),32,repmat(cmap(1,:),sum(sel1),1),'s','filled');
    %     s2 = scatter(newLatent{i}(logical(sel2),1),newLatent{i}(logical(sel2),2),24,repmat(cmap(2,:),sum(sel2),1),'d','filled');
    %     s3 = scatter(newLatent{i}(logical(sel3),1),newLatent{i}(logical(sel3),2),16,repmat(cmap(3,:),sum(sel3),1),'^','filled');
    %
    % Show prototype positions
    protOrg = output{i}.data{1}.prototypes(output{i}.data{1}.conceptSelection.id,:);
    [~, index]=ismember(protOrg,optima{i},'rows');
    protLat = newLatent{i}(index,:);
    s5 = scatter(protLat(1),protLat(2),256,'rd','filled');
    
    protOrg = output{i}.data{2}.prototypes(output{i}.data{2}.conceptSelection.id,:);
    [~, index]=ismember(protOrg,optima{i},'rows');
    protLat = newLatent{i}(index,:);
    s6 = scatter(protLat(1),protLat(2),256,cmap(1,:),'d','filled');
    
    legend([s1,s5,s4,s6,s2,s3],'Not Selected','Sel. Prototype 1','Iteration 1','Sel. Prototype 2','Iteration 2','Iteration 3','Location','SouthWest');
    ax = gca;
    ax.Visible = 'off';
end
%%
save_figures(fig, './', ['sail_vs_produqd_latent_'], 20, [6 5]);

%% Show examples
clear fig; set(0,'DefaultFigureWindowStyle','default');
cellStart = 8;
cellCent = 12;
cellEnd = 18;
selCells = [cellStart,cellStart;cellStart,cellCent;cellStart,cellEnd; ...
    cellCent,cellStart;cellCent,cellCent;cellCent,cellEnd; ...
    cellEnd,cellStart;cellEnd,cellCent;cellEnd,cellEnd];
numAlgs = 4;
cmap = parula(numAlgs);
lineStyles = {'-','--',':','-.'};
for i=1:size(selCells,1)
    fig(i) = figure(i);hold off;
    for alg=1:numAlgs
        geneExample = squeeze(output{alg}.data{end}.sail.predMap(end).genes(selCells(i,1),selCells(i,2),:));
        expressedExample = d.express(geneExample');
        h = fPlot(expressedExample,'Color',cmap(alg,:),'LineWidth', 4, 'LineStyle', lineStyles{alg});
        hold on;ax = gca;ax.Visible = 'off';
        axis equal; axis tight;
        if alg==numAlgs && i == size(selCells,1);legend('SAIL','PRODUQD1','PRODUQD2','PRODUQD3');end
    end
end
save_figures(fig, './', ['sail_vs_produqd_ontop_'], 20, [8 3]);

%% Show evolution of similarity space
cmap = parula(4);
for i=2:numPRODUQDRuns
    optima{i} = [];
    for j=1:p.numIterations
        optima{i} = [optima{i}; output{i}.data{j}.optima];
        numSolutions(i,j) = size(output{i}.data{j}.optima,1);        
        [~, latSeq{i,j}] = dimReducedClustering( optima{i}, 'tSNE', 2);
        if j > 1
            % Rotate latent space to improve visually matching it with first iteration
            A = latSeq{i,1};
            B = latSeq{i,j};
            centroid_A = mean(A);
            centroid_B = mean(B(1:size(latSeq{i,1},1),:));
            N = size(A,1);
            H = (B(1:size(latSeq{i,1},1),:) - repmat(centroid_B, N, 1))' * (A - repmat(centroid_A, N, 1));
            [U,S,V] = svd(H);
            R = V*U';
            t = -R*centroid_A';% + centroid_B';
            B2 = (R*B') + repmat(t, 1, size(B,1));
            latSeq{i,j} = B2';
        end            
    end
end

%% Visualization
clear fig;
cmap = parula(3);
for i=3:3%numPRODUQDRuns
    for j=1:p.numIterations
        fig(1) = figure(1);subplot(1,3,j);hold off;
        s1 = scatter(latSeq{i,j}(1:size(latSeq{i,1},1),1),latSeq{i,j}(1:size(latSeq{i,1},1),2),64,[0.7 0.7 0.7],'o','filled');hold on;
        s2 = scatter(latSeq{i,j}(size(latSeq{i,1},1)+1:end,1),latSeq{i,j}(size(latSeq{i,1},1)+1:end,2),64,cmap(2,:),'o','filled');hold on;
        if j == 3;s5 = scatter(latSeq{i,j}(sum(numSolutions(i,1:2))+1:end,1),latSeq{i,j}(sum(numSolutions(i,1:2))+1:end,2),64,cmap(3,:),'o','filled');hold on;end
        
        labels = output{i}.data{1}.conceptLabels;
        sel = logical(zeros(size(latSeq{i,1},1),1));
        firstClassIDs = [labels==output{i}.data{1}.conceptSelection.id];
        sel(1:length(firstClassIDs)) = firstClassIDs;
        s3 = scatter(latSeq{i,j}(sel,1),latSeq{i,j}(sel,2),64,cmap(1,:),'o','filled');hold on;
        
        protOrg = output{i}.data{1}.prototypes(output{i}.data{1}.conceptSelection.id,:);
        [~, index]=ismember(protOrg,optima{i},'rows');
        protLat = latSeq{i,j}(index,:);
        s4 = scatter(protLat(1),protLat(2),128,'ro','filled');
    
        axis([-50 50 -50 50]);
        legend([s1 s3 s4 s2], 'Not Selected', 'Selected Class', 'Prototype', 'New Optima 1', 'New Optima 2');
        if j == 3;legend([s1 s3 s4 s2 s5], 'Not Selected', 'Selected Class', 'Prototype', 'New Optima 1', 'New Optima 2');end
        ax = gca;ax.XTick = [];ax.YTick = [];
    end
end

save_figures(fig, './', ['sail_vs_produqd_evo_'], 20, [24 5]);
