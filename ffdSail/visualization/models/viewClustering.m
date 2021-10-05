clear;clc;

%domainname = 'PARSEC';systemInit;
%[experimentNames,data] = read_experiments('/scratch/ahagg2s/acq_parsec_aligned');
%[experimentNames,data] = read_experiments('/scratch/ahagg2s/acq_parsec_nonaligned');

domainname = 'FOILFFD';systemInit;
shownRuns = [5];
[experimentNames,data] = read_experiments('/scratch/ahagg2s/acq_ffd_1','runsToShow', shownRuns);


%% Extract statistics, dimensionality reduction and clustering
d = data{1}.d;p = data{1}.p;
d.alignedMap = false; %TODO get rid of this, should be saved in the experiment data anyway

disp(['Retrieving experiment statistics']);
[stats,p,d,collect] = getModelComparisonStats(data, p, d, experimentNames);
optima = stats.agg.optimaLocationsAndSamples;

disp(['Get optimal regions with t-SNE/DBSCAN']);
[estimatedLabels, reducedParSpace, ~, epsilon, coreneighbours, stats.clusters] = tsne_dbscan( optima, 'max_iter', 1000,  'perplexity', 200, 'coreneighbours', 3 );
disp([int2str(max(estimatedLabels))  ' clusters found, ' int2str(sum(estimatedLabels==0)) '/' int2str(size(reducedParSpace,1))  ' samples unassigned']);

disp(['Analyze optimal regions']);
[stats.clusters.analysis] = getClusterAnalysis(stats,p,d,collect, data, estimatedLabels, reducedParSpace);

%% Visualization: configure
set(0,'DefaultFigureWindowStyle','default')
clear fig;
dotsize = 16;

% Sort clusters by size
tbl = tabulate(estimatedLabels);
[clusterIDcounts,sortedLabelIDs] = sort(tbl(:,2),'descend');
sortedLabels = tbl(sortedLabelIDs,1);
sortedLabels(sortedLabels==0) = [];
clusterIDcounts(sortedLabels==0) = [];
% Select largest clusters
clusterSelection = 1:10;
showClusters = sortedLabels(clusterSelection);

% Created colormap, sorted by cluster size (first color is largest cluster)
colorsSelectedClusters = colorcube(length(showClusters));
colorsClusters = repmat([0.8 0.8 0.8],size(sortedLabels,1),1);
colorsClusters(showClusters,:) = colorsSelectedClusters;
%colorsClustersAllSamples = colorsClusters(clusterSelection,:);
%colorsFeatures = jet(625);colorsFeatures(1,:) = [0.8 0.8 0.8];
rotMat = reshape(1:625,25,25)';
%colorsFeaturesRot = colorsFeatures(rotMat(:),:);

fitnessBasedSizes = ceil(-[stats.clusters.analysis.sortedFitness(:); stats.clusters.analysis.sampleFitnesses']);
fitnessBasedSizes = fitnessBasedSizes(~stats.clusters.analysis.allNans);
fitnessBasedSizes = ceil(fitnessBasedSizes);
fitnessBasedSizes = 2.^fitnessBasedSizes;

%% Show clusters in feature map
for i=1:length(data)
    clusterIDMap = nan(625,1);
    clusterIDMap(~stats.agg.nanSolutions{i}) = estimatedLabels(stats.agg.optimaIDs{i});
    % Visualize
    clear fig;clf
    
    fig(1) = figure(1);hold off;
    
    [f,~,h] = viewMap(reshape(clusterIDMap,25,25),d);
    colormap(colorsClusters);
    caxis([1 max(estimatedLabels)]);

    h.Label.String = 'Largest clusters (color)';
    title(['Optimal regions in feature space']);
    %
    fig(2) = figure(2);hold off;
    % Show deselected samples
    h(1) = scatter(reducedParSpace(max(stats.agg.optimaIDs{end})+1:size(reducedParSpace,1),1),...
        reducedParSpace(max(stats.agg.optimaIDs{end})+1:size(reducedParSpace,1),2),...
        fitnessBasedSizes(max(stats.agg.optimaIDs{end})+1:size(reducedParSpace,1)),...
        [0.5 0.5 0.5],...
        'filled', 'MarkerFaceAlpha', 1);
    hold on;
    
    for m=1:length(stats.cfg.xpID)
        [~,selectSample] = ismember(estimatedLabels,showClusters,'rows');
        optima = stats.agg.optimaIDs{m};
        labels = selectSample(1:length(optima));
        labels(selectSample(1:length(optima))==0) = [];
        theseColors = colorsSelectedClusters(labels,:);
        h(2) = scatter(reducedParSpace(selectSample(1:length(optima))~=0,1),reducedParSpace(selectSample(1:length(optima))~=0,2),...
            fitnessBasedSizes(selectSample(1:length(optima))~=0),...
            theseColors,...
            's','filled',  'MarkerEdgeColor',[0 0 0], 'LineWidth',1);
    end
    ax = gca;
    ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; axis([min(reducedParSpace(:,1)) max(reducedParSpace(:,1)) min(reducedParSpace(:,2)) max(reducedParSpace(:,2))]);
    save_figures(fig, './', ['featToParMapping_'], 8, [3 3]);
end




%% Overlapping shapes of every cluster
close all;clf;
numPCExplained = 5; % Over how many principal components calculate "% of variance explained"
optima = stats.agg.optimaLocationsAndSamples;



axisRanges = [0 1 -0.4 0.4];
for ii=1:length(showClusters)
    % Get all nD samples in this cluster
    samples = optima(estimatedLabels==showClusters(ii),:);
    for ss=1:size(samples,1)
        sampleExpressed{ii}(ss,:,:) = d.express(samples(ss,:));
    end
end

for ii=1:length(showClusters)
    % Get all nD samples in this cluster
    samples = optima(estimatedLabels==showClusters(ii),:);
    size(samples,1)
end

% % PCA again
for ii=1:length(showClusters)
    % Get all nD samples in this cluster
    samples = optima(estimatedLabels==showClusters(ii),:);
    % Remove mean shape
    clSamples = samples - mean(samples);
    % PCA
    [coeff,score,latent,tsquared,explained,mu] = pca(clSamples);
    
    features = coeff(:,1)'*clSamples'; % Get first principal component
    [~,medFeatureID] = min(abs(features - median(features)));
    medFeature(ii,:) = features(medFeatureID);
    prototype{1}(ii,:,:) = d.express(samples(medFeatureID,:));
    
    features = coeff(:,1:2)'*clSamples'; % Get first principal component
    [~,medFeatureID] = min(sum(abs(features - median(features')')));
    medFeature(ii,:) = features(medFeatureID);
    prototype{2}(ii,:,:) = d.express(samples(medFeatureID,:));
    
end

% Get medoids of clusters
for ii=1:length(showClusters)
    % Get all nD samples in this cluster
    samples = optima(estimatedLabels==showClusters(ii),:);
    [~,c] = kmedoids(samples,1);
    prototype{3}(ii,:,:) = d.express(c);
end

%
for ii=1:length(showClusters)
    % Get all nD samples in this cluster
    samples = optima(estimatedLabels==showClusters(ii),:);
    fig(ii) = figure(ii);clf;hold on;
    for ss=1:size(samples,1)
        p1 = plot(squeeze(sampleExpressed{ii}(ss,1,:)),squeeze(sampleExpressed{ii}(ss,2,:)),'k','LineWidth',2);
        p1.Color(4) = 5/size(samples,1);
    end
    p3 = plot(squeeze(prototype{1}(ii,1,:)),squeeze(prototype{1}(ii,2,:)),'b--','LineWidth',4);
    hold on;
    p4 = plot(squeeze(prototype{2}(ii,1,:)),squeeze(prototype{2}(ii,2,:)),'b:','LineWidth',4);
    p5 = plot(squeeze(prototype{3}(ii,1,:)),squeeze(prototype{3}(ii,2,:)),'g-','LineWidth',4);
    axis equal;axis(axisRanges);grid on;
    title([int2str(size(samples,1)) ' Samples']);
    legend([p3 p4 p5],'Median 1st PC', 'Median 1st & 2nd PC', 'Cluster Medoid');
    drawnow;
end
disp(['% explained by first ' int2str(numPCExplained) ' principal components: ' ]);
disp(string(sum(explained(1:numPCExplained,:))));

%%
for i=1:length(showClusters)
    saveas(fig(i),['areTheyRepresentatives' int2str(i) '.png'])
end
%save_figures(fig, './', ['areTheyRepresentatives'], 12, [8 6]);
