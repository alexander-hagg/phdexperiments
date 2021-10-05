clear;clc;systemInit;
%%
allresults = {};output = {};
warning('off','all');
resultpath = '/scratch/ahagg2s/acq_gps_long4/1';
fsNodes = dir(resultpath);fsNodes([fsNodes.isdir]) = [];
% Load results from mat files
for f=1:length(fsNodes)
    if strcmp(fsNodes(f).name(end-2:end),'mat')
        output = load([fsNodes(f).folder '/' fsNodes(f).name]);
        allresults{end+1} = output.output;
    end
end
warning('on','all');

%% Load Sobol Set of same length
numSamples = size(allresults{1}.model{1}.trainInput,1);
if ~exist(['ffd_' int2str(numSamples) '.mat'],'file')
    sampleFileNames = ffd_GetSobolDatasets(d, 'rndSampleSizes', [numSamples]);
end
sobolSampleSet = load(['ffd_' int2str(numSamples) '.mat']);
value = [sobolSampleSet.cD_true' sobolSampleSet.cL_true'];
for iModel = 1:2 %cD and cL
    disp(['Train model ' int2str(iModel)]);
    sobolSampleSet.modelPred{iModel} = feval(['train' d.paramsPred{iModel}.name], sobolSampleSet.observation, value(:,iModel), d.paramsPred{iModel});
end
[predMap, predPercImprove] = createPredictionMap(sobolSampleSet.modelPred,p,d,'featureRes',[25 25]);
sobolSampleSet.predMap(numSamples) = predMap;
sobolSampleSet.predPercImprove(numSamples,:) = predPercImprove;
output = sobolSampleSet;
save('sobolSampleSet.mat','output', 'p', 'd', '-v7.3');
sobolSampleSet = getGroundTruth('sobolSampleSet.mat', 'domains/foilFFD', 'numSamples', 1000);
clear output;
%% Statistics
% Determine number of maps to evaluate
stats.cfg.prcNames = [25,75];
xpNames = {'GP_run_1','GP_run_2','GP_run_3', 'GP_run_4'};%, 'Sobol'};
xpNames = strrep(xpNames,'_','\_');
stats.cfg.runsToShow= 1:length(xpNames);
stats.cfg.binsize = 50;
stats.cfg.selectedReplicates = 1;
stats.cfg.overrideShownPeriod = 11;
d = allresults{stats.cfg.runsToShow(1)}.d;
p = allresults{stats.cfg.runsToShow(1)}.p;
[stats,p,d] = getModelComparisonStats(allresults, p, d, stats);

%%
groupSelected = [1:1:size(stats.analysis.medianPercentageMap,2)];
if groupSelected(1) ~=1; groupSelected = [1 groupSelected];end
if groupSelected(end) ~=size(stats.analysis.medianPercentageMap,2); groupSelected = [groupSelected size(stats.analysis.medianPercentageMap,2)];end
subgroups = ones(1,size(stats.analysis.medianPercentageMap,1));
clear fig;
fig(1) = figure(1);clf;hold off;

hndl = viewMedianPercentiles(100-100*stats.analysis.medianPercentageMap, 100-100*stats.analysis.prctilePercentageMap, ...
    'groupNames', stats.cfg.iItr, ...
    'subgroups', subgroups, ...
    'prcNames', string(stats.cfg.prcNames), ...
    'colorScheme', 'parula', ...
    'tickLen', 2 , ...
    'tickRelDis', 0.02, ...
    'groupDist', 4, ...
    'groupSelected', groupSelected);

xlabel('Number of Samples');
ylabel(['Fitness (% of max) (+ ' int2str(stats.cfg.prcNames(1)) '/' int2str(stats.cfg.prcNames(2)) ' perc.)']);
ax = gca;ax.XAxis.Scale = 'log';ax.YAxis.Scale = 'log';ax.YDir = 'reverse';
axis([0.9*p.nInitialSamples 1.1*p.nTotalSamples 0 100]);
ax.YTick = [1 10 50 100];
ax.YTickLabel = string((100-ax.YTick));

title(['Predicted Fitness, % of Max Fitness in Corresponding Bin']);
l = legend(hndl, xpNames, 'Location', 'NorthWest');


%%
fig(2) = figure(2);
hold off;

hndl = viewMedianPercentiles(100 - stats.analysis.output.MED, 100 - stats.analysis.output.PRC, ...
    'subgroups', ones(1,size(stats.analysis.medianPercentageMap,1)), ...
    'groupNames', stats.cfg.groups, ...
    'colorScheme', 'parula');
xlabel('Acquired samples (chronological order)');
ylabel(['Fitness (% of max) (+ ' int2str(stats.cfg.prcNames(1)) '/' int2str(stats.cfg.prcNames(2)) ' perc.)']);
legend(hndl,xpNames{stats.cfg.runsToShow}, 'Location', 'NorthWest');
ax = gca;ax.XAxis.Scale = 'log';ax.YAxis.Scale = 'log';ax.YDir = 'reverse';
axis([0.9*p.nInitialSamples 1.1*p.nTotalSamples 0 100]);
ax.YTick = [1 10 50 100];
ax.YTickLabel = string((100-ax.YTick));
title(['Acquired Samples, % of Max Fitness in Corresponding Bin']);
ax = gca;
ax.XAxis.Scale = 'log';

%%
fig(3) = figure(3);clf;hold off;
hndl = viewMedianPercentiles(stats.analysis.medianConfidenceMap, stats.analysis.prctileConfidenceMap, ...
    'groupNames', stats.cfg.iItr(1:end), ...
    'subgroups', subgroups, ...
    'prcNames', string(stats.cfg.prcNames), ...
    'colorScheme', 'parula', ...
    'tickLen', 2 , ...
    'tickRelDis', 0.02, ...
    'groupDist', 4, ...
    'groupSelected', groupSelected(1:end-1));

xlabel('Number of Samples');
ylabel(['Model Variance (+ ' int2str(stats.cfg.prcNames(1)) '/' int2str(stats.cfg.prcNames(2)) ' perc.)']);
title(['Median Bin Confidence']);
legend(hndl, xpNames, 'Location', 'NorthWest');


save_figures(fig, './', 'AcquisitionModelComparison', 14, [16 6]);

%%

%dimReduxMethods = {'PCA', 'LDA', 'MDS', 'ProbPCA', 'KernelPCA', 'GDA', 'tSNE', 'ManifoldChart'};
dimReduxMethods = {'PCA'};
for i=1:length(dimReduxMethods)
    samplesParams = [stats.analysis.in{1};sobolSampleSet.observation];
    ttt = compute_mapping(samplesParams, dimReduxMethods{i}, 2, 10, 1000);
    Y(i,:,:) = ttt;
end

Y = 10 * Y / max(abs(min(Y(:))),max(Y(:)));


clrsmap = [1:5:25];
clrsmap = repelem(clrsmap,1,5);
binColorV = repmat(clrsmap,25,1)/25;
binColorU = binColorV';
binColorW = ones(size(binColorV));
binColors = cat(3,binColorU,binColorV,binColorW);

fig(1) = figure(100);hold off;im = imshow(binColors);
title('Bins in Feature Space');

%% Add Sobol Set into the mix
stats.analysis.mapLinIndx{end+1} = sobolSampleSet.mapLinIndx;
xpNames{end+1} = 'Sobol';

%%
clear fig;
clrs = [1 0.1 0.3; 0.3 1 0.1; 0.1, 0.3, 1; 0.1, 0.1, 0.3];
set(0,'DefaultFigureWindowStyle','docked')
nlevels = 64;

xlin = -10:0.1:10; ylin = xlin;
[xq,yq] = meshgrid(xlin, ylin);

for i=1:length(xpNames)
    pickedSamples{i}= 1+p.nTotalSamples*(i-1):p.nTotalSamples*i;
end

clear density xi;
for i=1:length(xpNames)
    YY = squeeze(Y(1,pickedSamples{i},:));
    
    % Get Fitness Data
    if i < length(xpNames)
        fitnessSamples = -stats.analysis.trainingRealFit{1}(1+p.nTotalSamples*(i-1):p.nTotalSamples*i);
        maxFitnessSamples = stats.analysis.maxFitMap(stats.analysis.trainingBinIDs{1}(1+p.nTotalSamples*(i-1):p.nTotalSamples*i));
        fitnessSamples = fitnessSamples./maxFitnessSamples;
        fitq(i,:,:) = griddata(YY(:,1),YY(:,2),fitnessSamples,xq,yq, 'natural');
    else
        fitnessSamples = -sobolSampleSet.fitness_true;
        maxFitnessSamples = stats.analysis.maxFitMap(sobolSampleSet.mapLinIndx);
        fitnessSamples = fitnessSamples./maxFitnessSamples;
        fitq(i,:,:) = griddata(YY(:,1),YY(:,2),fitnessSamples,xq,yq, 'natural');
    end
    
    % Sampling Density
    % Estimate a continuous pdf from the discrete data
    [density(i,:) xi(i,:,:)]= ksdensity(YY,[xq(:) yq(:)], 'Bandwidth', 0.5);
    
    % Bin IDs
    binU(i,:,:) = griddata(YY(:,1),YY(:,2),double(binColorU(stats.analysis.mapLinIndx{i}(1:p.nTotalSamples))),xq,yq, 'nearest');
    binV(i,:,:) = griddata(YY(:,1),YY(:,2),double(binColorV(stats.analysis.mapLinIndx{i}(1:p.nTotalSamples))),xq,yq, 'nearest');
    binW(i,:,:) = griddata(YY(:,1),YY(:,2),double(binColorW(stats.analysis.mapLinIndx{i}(1:p.nTotalSamples))),xq,yq, 'nearest');
end
maxDensity = max(density(:));
density = density / maxDensity;

%%
clear fig;

for i=1:3%length(xpNames)
    thisDensity = reshape(density(i,:),sqrt(length(density)),sqrt(length(density(i,:))));
    minDensity = (thisDensity>0.1);    
    
    fig((i-1)*3 + 1) = figure((i-1)*3 + 1);hold off;
    [~,h] = contourf(xq,yq,minDensity.*squeeze(fitq(i,:,:)),nlevels);
    set(h,'LineColor','none');
    c = colorbar;c.Label.String = 'Fitness';
    caxis([0 1]);axis equal;axis([min(xlin) max(xlin) min(ylin) max(ylin)]);set(gca,'XTick',[]); set(gca,'YTick',[]);
    ylabel(['Acquisition: ' strrep(xpNames{i},':GP','')]);
    title('Rel. Bin Fitness');
    
    fig((i-1)*3 + 2) = figure((i-1)*3 + 2);hold off;
    [~,h] = contourf(xq,yq,minDensity.*thisDensity,nlevels);
    set(h,'LineColor','none')
    c = colorbar;c.Label.String = 'Sampling Prob. Density';
    caxis([0 1]);axis equal;axis([min(xlin) max(xlin) min(ylin) max(ylin)]);set(gca,'XTick',[]); set(gca,'YTick',[]);
    title('Sampling Density (pdf)');
    
    fig((i-1)*3 + 3) = figure((i-1)*3 + 3);hold off;
    imagesc(cat(3,flipud(squeeze(binU(i,:,:))),flipud(squeeze(binV(i,:,:))),flipud(squeeze(binW(i,:,:)))));
    set(h,'LineColor','none');
    axis equal;axis tight;
    title('Bin IDs');
    drawnow;
end
%save_figures(fig, './', ['clusters'], 14, [6 6]);

%%
    fig(1) = figure(1);hold off;
    [~,h] = contourf(xq,yq,squeeze(nanstd(fitq(1:4,:,:))),nlevels);
    set(h,'LineColor','none');
    c = colorbar;c.Label.String = 'Fitness';
    caxis([0 0.3]);axis equal;axis([min(xlin) max(xlin) min(ylin) max(ylin)]);set(gca,'XTick',[]); set(gca,'YTick',[]);
    ylabel(['Acquisition: ' strrep(xpNames{i},':GP','')]);
    title('\sigma Rel. Bin Fitness');
    
    %
    fig(2) = figure(2);hold off;
    thisDensity = reshape(std(density(1:4,:)),sqrt(length(density)),sqrt(length(density(1,:))));
    [~,h] = contourf(xq,yq,thisDensity,nlevels);
    set(h,'LineColor','none')
    c = colorbar;c.Label.String = 'Sampling Prob. Density';
    caxis([0 0.3]);axis equal;axis([min(xlin) max(xlin) min(ylin) max(ylin)]);set(gca,'XTick',[]); set(gca,'YTick',[]);
    title('\sigma Sampling Density');
    
    ggg = griddata(YY(:,1),YY(:,2),double(stats.analysis.mapLinIndx{i}(1:p.nTotalSamples)),xq,yq, 'nearest');    
    fig(3) = figure(3);hold off;
    [~,h] = contourf(xq,yq,ggg,nlevels);
    %imagesc(cat(3,flipud(squeeze(std(binU))),flipud(squeeze(std(binV))),flipud(squeeze(std(binW)))));
    set(h,'LineColor','none');
    axis equal;axis tight;
    title('\sigma Bin IDs');
    drawnow;

%%
%clear fig;
fig(end+1) = figure(99);hold off;
clrs = [1 0.1 0.3; 0.3 1 0.1; 0.1, 0.3, 1; 0.1, 0.1, 0.3];
clrs = repelem(clrs,p.nTotalSamples,1);
for j=1:length(dimReduxMethods)
    subplot(ceil(sqrt(length(dimReduxMethods))),floor(sqrt(length(dimReduxMethods))),j);
    for repID=stats.cfg.selectedReplicates
        rep = stats.cfg.selectedReplicates(repID);
        for i=1:length(xpNames)
            fitnessSamples = -stats.analysis.trainingRealFit{rep}(1+p.nTotalSamples*(i-1):p.nTotalSamples*i);
            maxFitnessSamples = stats.analysis.maxFitMap(stats.analysis.trainingBinIDs{rep}(1+p.nTotalSamples*(i-1):p.nTotalSamples*i));
            fitnessSamples = fitnessSamples./maxFitnessSamples;
            siz = fitnessSamples*32;
            h = scatter(squeeze( Y(j,1+p.nTotalSamples*(i-1)+2*p.nTotalSamples*(rep-1):p.nTotalSamples*i+2*p.nTotalSamples*(rep-1),1) ), ...
                squeeze( Y(j,1+p.nTotalSamples*(i-1)+2*p.nTotalSamples*(rep-1):p.nTotalSamples*i+2*p.nTotalSamples*(rep-1),2) ), ...
                32, clrs(1+p.nTotalSamples*(i-1):p.nTotalSamples*i,:), ...
                'filled', ...
                'MarkerFaceAlpha', .7, ...
                'MarkerEdgeAlpha', .7 );
            set(h, 'SizeData', siz)
            hold on;
        end
    end
    axis([min(Y(j,:,1)) max(Y(j,:,1)) min(Y(j,:,2)) max(Y(j,:,2))]);
    l = legend(xpNames{stats.cfg.runsToShow}, 'location', 'southeast');
    title(['Distribution of samples in parameter space (' dimReduxMethods{j} ')' ]);
    hold on;
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    axis equal;
    drawnow;
end

%save_figures(fig, './', ['distribution_samples2'], 14, [10 10]);