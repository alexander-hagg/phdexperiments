clear;clc;
numSelConcepts = 5;
domainname = 'FOILFFD';systemInit; xpfoldername = '/scratch/ahagg2s/acq_ffd_1'; shownRuns = [1:30];

%% Run SAIL 5 times selecting largest clusters based on first 500 PE
for run=1:length(shownRuns)
    % Read 1000PE SAIL runs
    [experimentNames,data] = read_experiments(xpfoldername,'runsToShow', shownRuns(run));
    d = data{1}.d;d.tmpdir = '/tmp/xfoil';d.fitnessPenalty_Area = true;
    p = data{1}.p;d.alignedMap = false; %TODO get rid of this, should be saved in the experiment data anyway
    
    % Edit Hyperparameters
    p.nInitialSamples   = 500;
    p.nTotalSamples     = 1000;
    
    % Prepare samples
    observation = data{1}.samples.genes(1:p.nInitialSamples,:);
    value = [data{1}.samples.cD_true(1:p.nInitialSamples,:),data{1}.samples.cL_true(1:p.nInitialSamples,:)];
    d.initialSampleSource = ['run_' int2str(run) '_data'];
    d.loadInitialSamples = true;
    save(d.initialSampleSource,'observation','value');
    
    % Prepare optima found after 500 PE
    optima = reshape(data{1}.predMap(end/2).genes,625,10);
    optima(any(isnan(optima')),:) = [];
    
    % Extract clusters
    [o{run}.estimatedLabels, o{run}.reducedParSpace, ~, o{run}.epsilon, o{run}.coreneighbours, o{run}.clusters] = ...
        dimReducedClustering( optima, 'tSNE', 2);
    
    [~,sortedClusters] = sort(o{run}.clusters.sizes(2:end), 'descend'); % Exclude nonassigned clusters (first entry, ID = 0)
    
    clusters = sortedClusters(1:numSelConcepts);
    
    % Run SAIL
    for c = 1:length(clusters)
        % Limit SAIL to cluster
        % Cluster ID
        p.concept.id = clusters(c);
        options = statset('Display','final');
        p.concept.members = optima(o{run}.estimatedLabels==p.concept.id,:);
        % Multivariate Gaussian Distribution determines mutation bias (soft
        % constraints
        %% TODO: hard constraints
        p.concept.value = gmdistribution.fit(p.concept.members,1,'Options',options);
        
        runTime = tic;
        o{run}.sail{c} = sail(p,d);
        disp(['Runtime: ' seconds2human(toc(runTime))]);
        disp('Get ground truth of prediction maps');
        runTimeGT = tic;
        disp('Get ground truth of prediction maps');
        o{run}.sail{c}.predMap = getGroundTruthMapArray(o{run}.sail{c}.predMap, d, 'dummy', false);
        disp(['Runtime: ' seconds2human(toc(runTime)) ', of which ground truth took: ' seconds2human(toc(runTimeGT))]);
        
    end
    
end



%% Check whether new samples are in the cluster
postMortem_newOptima = reshape(o{run}.sail{c}.predMap(end).genes,625,10);
postMortem_newOptima(any(isnan(postMortem_newOptima')),:) = [];
[postMortem_labels, postMortem_redX, ~, ~, ~, postMortem_clusters] = dimReducedClustering( [optima;postMortem_newOptima], 'tSNE', 2);
postMortem_oldLabels = o{run}.estimatedLabels;
postMortem_oldRedX = o{run}.reducedParSpace;
postMortem_oldClusters = o{run}.clusters;

c = 4;
newLabelsOldCluster = unique(postMortem_labels(o{run}.estimatedLabels==clusters(c)));
newLabelsOldCluster(newLabelsOldCluster==0) = [];

oldClusterIndivs = o{run}.estimatedLabels==clusters(c);
newClustersContainOld = ismember(postMortem_labels,newLabelsOldCluster);

%
set(0,'DefaultFigureWindowStyle','docked')

clear fig;
fig(1) = figure(1);hold off;
scatter(o{run}.reducedParSpace(:,1), o{run}.reducedParSpace(:,2),32,[1 0 0],'filled');
hold on;
scatter(o{run}.reducedParSpace(o{run}.estimatedLabels==clusters(c),1), o{run}.reducedParSpace(o{run}.estimatedLabels==clusters(c),2),32,[0 0 1],'filled');
legend('Old Optima', 'Selected Concept','Location','SouthEast');
title('Original tSNE dimensions');
ax = gca;ax.XAxis.Visible = 'off';ax.YAxis.Visible = 'off';

for i=2:3
    fig(i) = figure(i);hold off;
    scatter(postMortem_redX(1:size(o{run}.estimatedLabels),1),postMortem_redX(1:size(o{run}.estimatedLabels),2),32,[1 0 0],'filled');
    hold on;
    %if i>1; scatter(postMortem_redX(newClustersContainOld,1),postMortem_redX(newClustersContainOld,2),64,[0 0 1],'filled');end;
    if i>1; scatter(postMortem_redX(oldClusterIndivs,1),postMortem_redX(oldClusterIndivs,2),64,[0 0 1],'filled');end;
    if i>2;scatter(postMortem_redX(size(o{run}.estimatedLabels)+1:end,1),postMortem_redX(size(o{run}.estimatedLabels)+1:end,2),32,[0 1 0],'filled');end;
    legend('Old Optima', 'Selected Concept', 'New Optima','Location','SouthEast');
    title('New tSNE dimensions');
    ax = gca;ax.XAxis.Visible = 'off';ax.YAxis.Visible = 'off';
end


fig(4) = figure(4);hold off;
scatter(postMortem_redX(1:size(o{run}.estimatedLabels),1),postMortem_redX(1:size(o{run}.estimatedLabels),2),32,[1 0 0],'filled');
hold on;
newClusterIDs = postMortem_labels;
newUniqueClusterIDs = unique(newClusterIDs);
newUniqueClusterIDs(newUniqueClusterIDs==0) = [];
colorsSelectedClusters = colorcube(length(newUniqueClusterIDs)+10);
colorsSelectedClusters = colorsSelectedClusters(1:end-10,:);
colorsClusters = repmat([0.8 0.8 0.8],size(newUniqueClusterIDs,1),1);
colorsClusters = colorsSelectedClusters;

newRedX = postMortem_redX(size(o{run}.estimatedLabels)+1:end,:);
newClusterIDs = newClusterIDs(size(o{run}.estimatedLabels)+1:end);
newClusterIDsNoZeros = newClusterIDs(newClusterIDs~=0);
scatter(newRedX(newClusterIDs~=0,1),newRedX(newClusterIDs~=0,2),32,colorsClusters(newClusterIDsNoZeros,:),'filled');
title('New Subconcepts');
ax = gca;ax.XAxis.Visible = 'off';ax.YAxis.Visible = 'off';
%save_figures(fig, './', ['prodigi_' int2str(run) '_cluster' int2str(c)], 16, [7 6]);


%% Show selected and new prototype(s)
%fig(5) = figure(5);hold off;
% Sort clusters by size
tbl = tabulate(newClusterIDs);[clusterIDcounts,sortedLabelIDs] = sort(tbl(:,2),'descend');
sortedLabels = tbl(sortedLabelIDs,1);sortedLabels(sortedLabels==0) = [];clusterIDcounts(sortedLabels==0) = [];
% Select largest clusters
clusterSelection = 1:5;

fig = viewPrototypes([optima;postMortem_newOptima], sortedLabels(clusterSelection), postMortem_labels, d);

save_figures(fig, './', ['prodigi_' int2str(run) '_cluster' int2str(c) '_prototypes'], 12, [7 6]);

oldClusterIDs = o{run}.estimatedLabels;
fig = viewPrototypes(optima, clusters([1:5]), oldClusterIDs, d);
save_figures(fig, './', ['prodigi_' int2str(run) '_cluster' int2str(c) '_SELECTEDprototype'], 12, [7 6]);
    
    
%% Visualization
set(0,'DefaultFigureWindowStyle','docked')

axisFIT = [-5 -3];
axisDIS = [0 2];
axisERR = [0 5];
% Show SAIL fitness comparison
for run=1:1%length(shownRuns)
    for c = 1
%    for c = 1:length(clusters)
        [idx,center] = kmedoids(o{run}.sail{c}.p.concept.members,1);
        
        fig(1) = figure(1);
        viewMap(data{1}.predMap(end/2).fitness_true,d);title(['SAIL ' int2str(p.nInitialSamples) 'PE']);
        caxis(axisFIT);
        fits(1,:) = reshape(data{1}.predMap(end/2).fitness_true,numel(data{1}.predMap(end/2).fitness_true),1);
        
        fig(2) = figure(2);
        viewMap(data{1}.predMap(end).fitness_true,d);title(['SAIL ' int2str(p.nTotalSamples) 'PE']);
        caxis(axisFIT);
        fits(2,:) = reshape(data{1}.predMap(end).fitness_true,numel(data{1}.predMap(end).fitness_true),1);
        
        fig(3) = figure(3);
        viewMap(o{run}.sail{c}.predMap(end).fitness_true,d);title(['PRODIGI ' int2str(p.nTotalSamples) 'PE']);
        caxis(axisFIT);
        fits(3,:) = reshape(o{run}.sail{c}.predMap(end).fitness_true,numel(o{run}.sail{c}.predMap(end).fitness_true),1);
        
        fig(4) = figure(4);
        dist1 = pdist2(center,reshape(data{1}.predMap(end/2).genes,625,10));
        viewMap(reshape(dist1,25,25),d);title('Distance to prototype');
        caxis(axisDIS);
        dists(1,:) = reshape(dist1,numel(dist1),1);

        fig(5) = figure(5);
        dist1 = pdist2(center,reshape(data{1}.predMap(end).genes,625,10));
        viewMap(reshape(dist1,25,25),d);title('Distance to prototype');
        caxis(axisDIS);
        dists(2,:) = reshape(dist1,numel(dist1),1);
        
        fig(6) = figure(6);
        dist1 = pdist2(center,reshape(o{run}.sail{c}.predMap(end).genes,625,10));
        viewMap(reshape(dist1,25,25),d);title('Distance to prototype');
        caxis(axisDIS);
        dists(3,:) = reshape(dist1,numel(dist1),1);
        
        
        fig(7) = figure(7);
        %err = (data{1}.predMap(end/2).fitness_true - data{1}.predMap(end/2).fitness).^2;
        err = abs(data{1}.predMap(end/2).fitness_true - data{1}.predMap(end/2).fitness);
        viewMap(reshape(err,25,25),d);title('Abs. Error');
        caxis(axisERR);
        errors(1,:) = reshape(err,numel(err),1);
        
        fig(8) = figure(8);
        %err = (data{1}.predMap(end).fitness_true - data{1}.predMap(end).fitness).^2;
        err = abs(data{1}.predMap(end).fitness_true - data{1}.predMap(end).fitness);
        viewMap(reshape(err,25,25),d);title('Abs. Error');
        caxis(axisERR);
        errors(2,:) = reshape(err,numel(err),1);
        
        fig(9) = figure(9);
        %err = (o{run}.sail{c}.predMap(end).fitness_true - o{run}.sail{c}.predMap(end).fitness).^2;
        err = abs(o{run}.sail{c}.predMap(end).fitness_true - o{run}.sail{c}.predMap(end).fitness);
        viewMap(reshape(err,25,25),d);title('Abs. Error');
        caxis(axisERR);
        errors(3,:) = reshape(err,numel(err),1);


        drawnow;
        %save_figures(fig, './', ['sail_on_cluster_run' int2str(run) '_cluster' int2str(c)], 16, [6 6]);
    end
end

%% Boxplots
clear fig;%close all;
set(0,'DefaultFigureWindowStyle','default')

fig(1) = figure(1);
boxplot(fits', 'PlotStyle','compact', 'Colors', 'k');grid on;axis([0.5 3.5 -5.5 -3]);
ax = gca;
ax.XTickLabels = {'SAIL500', 'SAIL1000', 'PRODIGI1000'};
ylabel('True Fitness');


fig(2) = figure(2);
boxplot(dists', 'PlotStyle','compact', 'Colors', 'k');grid on;axis([0.5 3.5 -0.01 2.5]);
ax = gca;
ax.XTickLabels = {'SAIL500', 'SAIL1000', 'PRODIGI1000'};
ylabel('Distance to Prototype');

fig(3) = figure(3);
boxplot(errors', 'PlotStyle','compact', 'Colors', 'k');grid on;axis([0.5 3.5 -0.1 0.5]);
ax = gca;
ax.XTickLabels = {'SAIL500', 'SAIL1000', 'PRODIGI1000'};
ylabel('Absolute Error');

save_figures(fig, './', ['sail_on_cluster_errors'], 18, [3 6]);


%% Show Model Accuracy
for run=1:1%length(shownRuns)
    for c = 1:length(clusters)
        drawnow;
        save_figures(fig, './', ['sail_on_cluster_run' int2str(run) '_cluster' int2str(c)], 12, [4 7]);
    end
end


%% ANIMATED
for run=1:1%length(shownRuns)
    for c = 1:length(clusters)
        [idx,center] = kmedoids(o{run}.sail{c}.p.concept.members,1);
        
        file1 = 'sail_clusterdistance.gif';
        file2 = 'prodigi_clusterdistance.gif';
        
        for i=10:20
            
            fig(1) = figure(1);
            dist1 = pdist2(center,reshape(data{1}.predMap(i).genes,625,10));
            viewMap(reshape(dist1,25,25),d);%title('Distance to cluster center');
            caxis([0 2]);
            title(['SAIL Cluster ' int2str(c) ', ' int2str(i*50) 'PE']);
            drawnow;
            frame = getframe(fig(1));
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i == 10
                imwrite(imind,cm,[file1 '_cluster_' int2str(c)],'gif', 'Loopcount', 0, 'DelayTime', 0.2);
            else
                imwrite(imind,cm,[file1 '_cluster_' int2str(c)],'gif','WriteMode','append');
            end
            
            fig(2) = figure(2);
            dist1 = pdist2(center,reshape(o{run}.sail{c}.predMap(i*50).genes,625,10));
            viewMap(reshape(dist1,25,25),d);%title('Distance to cluster center');
            caxis([0 2]);
            title(['PRODIGI Cluster ' int2str(c) ', ' int2str(i*50) 'PE']);
            
            drawnow;
            frame = getframe(fig(2));
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if i == 10
                imwrite(imind,cm,[file2 '_cluster_' int2str(c)],'gif', 'Loopcount', 0, 'DelayTime', 0.2);
            else
                imwrite(imind,cm,[file2 '_cluster_' int2str(c)],'gif','WriteMode','append');
            end
            
            
            %save_figures(fig, '.', ['DISTANCE_cluster' int2str(c) '_iter_' int2str(i*50) '_'], 12, [3 3]);
            
        end
    end
    
end