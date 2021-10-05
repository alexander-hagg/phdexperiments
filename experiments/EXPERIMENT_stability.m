clear;clc;

numShownClusters = 10;
numDims = NaN; % Set to NaN to use intrinsic dimensionality
dimReduxMethods = {'none', 'PCA', 'KernelPCA', 'Isomap', 'tSNE', 'Autoencoder'};
%domainname = 'PARSEC';systemInit; xpfoldername = '/scratch/ahagg2s/acq_parsec_aligned'; shownRuns = [1 3:12 14:16]; %domainname = 'PARSEC';systemInit; xpfoldername = '/scratch/ahagg2s/acq_parsec_nonaligned'; shownRuns = [2:12 14];
domainname = 'FOILFFD';systemInit; xpfoldername = '/scratch/ahagg2s/acq_ffd_1'; shownRuns = [1,1,1,1,1,1,1,1,1,1];


for iii=1:length(shownRuns)
    run = shownRuns(iii);
    [experimentNames,data] = read_experiments(xpfoldername,'runsToShow', run);
    
    %% Extract statistics, dimensionality reduction and clustering
    d = data{1}.d;p = data{1}.p;d.alignedMap = false; %TODO get rid of this, should be saved in the experiment data anyway
    disp(['Retrieving experiment statistics']);
    [stats{iii},p,d,collect] = getModelComparisonStats(data, p, d, experimentNames, 'includesamples', false);
    
    for i=1:length(dimReduxMethods)
        disp('===============================================');
        disp(['Get optimal regions with ' dimReduxMethods{i} '/DBSCAN']);
        disp('===============================================');
        runData{iii,i}.method = dimReduxMethods{i};
        if isnan(numDims) || ~strcmp(dimReduxMethods{i},'tSNE')
            runData{iii,i}.numDims = ceil(intrinsic_dim(stats{iii}.agg.optimaLocations, 'MLE'));
        else
            runData{iii,i}.numDims = 2;
        end
        
        [runData{iii,i}.estimatedLabels, runData{iii,i}.reducedParSpace, ~, runData{iii,i}.epsilon, runData{iii,i}.coreneighbours, runData{iii,i}.clusters] = ...
            dimReducedClustering( stats{iii}.agg.optimaLocations, dimReduxMethods{i}, runData{iii,i}.numDims);
        disp([int2str(max(runData{iii,i}.estimatedLabels))  ' clusters found, ' int2str(sum(runData{iii,i}.estimatedLabels==0)) '/' int2str(size(runData{iii,i}.reducedParSpace,1))  ' samples unassigned']);
        
        disp(['Analyze optimal regions']);
        [runData{iii,i}.clusters.analysis] = getClusterAnalysis(stats{iii},p,d,collect, data, runData{iii,i}.clusters, runData{iii,i}.estimatedLabels, runData{iii,i}.reducedParSpace);
        
    end
end

%% Analyze Quality of Dimensionality Reduction
for iii=1:length(shownRuns)
    run = shownRuns(iii);
    sOrg = stats{iii}.agg.optimaLocations;
    for i=1:length(dimReduxMethods)
        sRed = runData{iii,i}.reducedParSpace;
        %Dd = pdist2(sOrg,sOrg); Dd = Dd ./ max(Dd(:)); % Normalized distance matrix original space
        %Dq = pdist2(sRed,sRed); Dq = Dq ./ max(Dq(:)); % Normalized distance matrix reduced space
        
        %% Normalized Shepard-Kruskal scaling
        %eks(iii,i) = norm((Dd-Dq),'fro').^2;
        
        %% Quality of k-neighbourhood ranking
        %enx(iii,i) = m_enx(Dd,Dq)
        
        %% Silhouette Coefficient
        %silh(iii,i) = m_silhouette(stats{iii}.agg.optimaLocations,runData{iii,i}.estimatedLabels);
        
        %% Calinski-Harabasz index
        %calhar(iii,i) = m_calinski(stats{iii}.agg.optimaLocations,runData{iii,i}.estimatedLabels);
        
        %% G+ index
        Gplus(iii,i) = m_gplus(pdist2(sRed,sRed), runData{iii,i}.estimatedLabels)
        
        %% Number of clusters
        numClusters(iii,i) = length(unique(runData{iii,i}.estimatedLabels));
    end
end
% save_figures(fig, './', ['coranking_matrix'], 12, [6 5]);  
%median(eks)
%median(enx)
%median(silh)
%median(calhar)

%%
for iii=1:length(shownRuns)
    for i=1:length(dimReduxMethods)
        csizes(iii,i) = length(runData{iii,i}.clusters.sizes);
    end
end

figure(1);
subplot(2,1,1);hold off;
plot(csizes,'LineWidth',2);
hold on;
plot(mean(csizes'),':','LineWidth',4);
grid on; legend(dimReduxMethods);
subplot(2,1,2);hold off;

%%
for run=1:length(shownRuns)
    for i=1:length(dimReduxMethods)
        for c=0:max(runData{run,i}.estimatedLabels)
            cX = runData{run,i}.reducedParSpace(runData{run,i}.estimatedLabels==c,:);
            X = stats{run}.agg.optimaLocations(runData{run,i}.estimatedLabels==c,:);
            mstd(run,i) = mean(std(cX));
        end
    end
end

%% new concept metric
for run=1:max(shownRuns)
    for i=1:length(dimReduxMethods)
        
    end
end


%% Visualization
for run=1:max(shownRuns)
    
    for i=1:length(dimReduxMethods)
        close all;
        % Visualization: configure
        set(0,'DefaultFigureWindowStyle','default')
        clear fig;
        dotsize = 16;
        
        % Sort clusters by size
        tbl = tabulate(runData{run,i}.estimatedLabels);
        [clusterIDcounts,sortedLabelIDs] = sort(tbl(:,2),'descend');
        sortedLabels = tbl(sortedLabelIDs,1);
        sortedLabels(sortedLabels==0) = [];
        clusterIDcounts(sortedLabels==0) = [];
        % Select largest clusters
        clusterSelection = 1:min(numShownClusters,max(sortedLabels));
        runData{run,i}.showClusters = sortedLabels(clusterSelection);
        showClusters = runData{run,i}.showClusters;
        
        % Created colormap, sorted by cluster size (first color is largest cluster)
        runData{run,i}.colorsSelectedClusters = colorcube(length(showClusters)+10);
        runData{run,i}.colorsSelectedClusters = runData{run,i}.colorsSelectedClusters(1:end-10,:);
        runData{run,i}.colorsClusters = repmat([0.8 0.8 0.8],size(sortedLabels,1),1);
        runData{run,i}.colorsClusters(showClusters,:) = runData{run,i}.colorsSelectedClusters;
        
        fitnessBasedSizes = ceil(-[runData{run,i}.clusters.analysis.sortedFitness(:); runData{run,i}.clusters.analysis.sampleFitnesses']);
        fitnessBasedSizes = fitnessBasedSizes(~runData{run,i}.clusters.analysis.allNans);
        fitnessBasedSizes = ceil(fitnessBasedSizes);
        runData{run,i}.fitnessBasedSizes = 2.^fitnessBasedSizes;
        
        showClusters = runData{run,i}.showClusters;
        estimatedLabels = runData{run,i}.estimatedLabels;
        
%         %% Parameter ranges
%         close all;clear fig;
%         optima = stats{run}.agg.optimaLocations;
%         fig(1) = figure(1);
%         for ii=1:length(showClusters)
%             % Get all nD samples in this cluster
%             samples = optima(estimatedLabels==showClusters(ii),:);
%             subplot(length(showClusters),1,ii);
%             boxplot(samples);
%             grid on;
%             ylabel(['Cluster ' int2str(ii)]);
%             axis([0.5 length(showClusters)+0.5 0 1]);
%         end
%         
%         save_figures(fig, './', ['parameterRanges_run' int2str(run) '_' dimReduxMethods{i}], 14, [16 12]);

        %% Show cluster prototypes
        fig = viewPrototypes(stats{run}.agg.optimaLocations, showClusters, estimatedLabels, d);
        %save_figures(fig, './', ['prototypes_run' int2str(run) '_' dimReduxMethods{i}], 12, [3 3], 300);
        for ff=1:length(fig)
            saveas(fig(ff),['prototypes_run' int2str(run) '_' dimReduxMethods{i} '_' int2str(ff) '.png'])
        end
        
        
        %% Show clusters in feature map
%         close all;clear fig;fig(1) = figure(1);
%         
%         clusterIDMap = nan(625,1);clusterIDMap(~stats{run}.agg.nanSolutions{1}) = estimatedLabels(stats{run}.agg.optimaIDs{1});
%         [f,~,h] = viewMap(reshape(clusterIDMap,25,25),d);
%         colormap(runData{run,i}.colorsClusters);
%         caxis([0.9 max(estimatedLabels)]);
%         h.Label.String = 'Largest clusters (color)';
%         title([dimReduxMethods{i} ': Optimal regions in feature space']);
%         
%         fig(2) = figure(2);
%         % Show deselected samples
%         reducedParSpace = runData{run,i}.reducedParSpace;
%         
%         h(1) = scatter(reducedParSpace(max(stats{run}.agg.optimaIDs{end})+1:size(reducedParSpace,1),1),...
%             reducedParSpace(max(stats{run}.agg.optimaIDs{end})+1:size(reducedParSpace,1),2),...
%             runData{run,i}.fitnessBasedSizes(max(stats{run}.agg.optimaIDs{end})+1:size(reducedParSpace,1)),...
%             [0.5 0.5 0.5],...
%             'filled', 'MarkerFaceAlpha', 1);
%         hold on;
%         for m=1:length(stats{run}.cfg.xpID)
%             [~,selectSample] = ismember(estimatedLabels,showClusters,'rows');
%             optima = stats{run}.agg.optimaIDs{m};
%             labels = selectSample(1:length(optima));
%             labels(selectSample(1:length(optima))==0) = [];
%             theseColors = runData{run,i}.colorsSelectedClusters(labels,:);
%             h(2) = scatter(reducedParSpace(selectSample(1:length(optima))~=0,1),reducedParSpace(selectSample(1:length(optima))~=0,2),...
%                 runData{run,i}.fitnessBasedSizes(selectSample(1:length(optima))~=0),...
%                 theseColors,...
%                 's','filled',  'MarkerEdgeColor',[0 0 0], 'LineWidth',1);
%         end
%         ax = gca;
%         ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off'; axis([min(reducedParSpace(:,1)) max(reducedParSpace(:,1)) min(reducedParSpace(:,2)) max(reducedParSpace(:,2))]);
%         title([dimReduxMethods{i} ': Clusters in ' int2str(numDims) 'D']);
%         save_figures(fig, './', ['featToParMapping_run' int2str(run) '_' dimReduxMethods{i}], 14, [6 6]);

    end
end
