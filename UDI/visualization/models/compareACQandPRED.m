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

%% Statistics
stats.cfg.prcNames = [25,75];
xpNames = strrep({'GP_run_1','GP_run_2','GP_run_3', 'GP_run_4'},'_','\_');
stats.cfg.runsToShow= 1:length(xpNames);
stats.cfg.binsize = 50;
stats.cfg.selectedReplicates = 1;
stats.cfg.overrideShownPeriod = 11;
d = allresults{stats.cfg.runsToShow(1)}.d;p = allresults{stats.cfg.runsToShow(1)}.p;
[stats,p,d] = getModelComparisonStats(allresults, p, d, stats);

allIndividuals = [];
allIndividuals = stats.analysis.in{1};
nSamples = size(allIndividuals,1)/length(xpNames);
acquisitionSelection = [0 nSamples 2*nSamples 3*nSamples 4*nSamples];

predictedIndividuals = stats.analysis.output.samples;
for i=1:length(xpNames)
    s = squeeze(predictedIndividuals(i,:,:,:));
    s = reshape(s,size(s,1)*size(s,2),size(s,3));
    
    allIndividuals = [allIndividuals; s(~isnan(stats.analysis.output.fitness_true(i,:)'),:)];
    acquisitionSelection(end+1) = acquisitionSelection(end) + sum(~isnan(stats.analysis.output.fitness_true(i,:)'));
end

stats.analysis.samples.reduced = compute_mapping(allIndividuals, 'PCA', 2, 10, 1000);

%%
nlevels = 32;
clear fitq density xi
xlin = -1:0.01:1; ylin = xlin;
[xq,yq] = meshgrid(xlin, ylin);
for j=1:2 % Acquired and Predicted Samples
for i=1:length(xpNames)
    YY = stats.analysis.samples.reduced((acquisitionSelection(i+(j-1)*length(xpNames))+1):acquisitionSelection((j-1)*length(xpNames)+i+1),:);
    
    % Get Fitness Data
    if j==1
        fitnessSamples = -stats.analysis.trainingRealFit{1}(1+p.nTotalSamples*(i-1):p.nTotalSamples*i);
        maxfit = stats.analysis.binsTraining{1}(1+p.nTotalSamples*(i-1):p.nTotalSamples*i);
        fitnessSamples = fitnessSamples./maxfit;
    else
        fit = -stats.analysis.output.fitness_true;
        fit = fit(i,:,:); fit = fit(:);
        fitnessSamples = fit(~isnan(stats.analysis.output.fitness_true(i,:)));
        fitnessSamples = fitnessSamples./stats.analysis.maxFitMap((~isnan(stats.analysis.output.fitness_true(i,:))))';    
    end
    
    fitq(j,i,:,:) = griddata(YY(:,1),YY(:,2),fitnessSamples,xq,yq, 'nearest');

    % Sampling Density
    % Estimate a continuous pdf from the discrete data
    [density(j,i,:) xi(j,i,:,:)]= ksdensity(YY,[xq(:) yq(:)], 'Bandwidth', 0.1);
    
end
end

%%


%%
% clear fig;
% 
% for i=1:length(xpNames)
%     minDensityThreshold = 0.1;
%     % Fitness and density acquired samples
%     thisDensity = reshape(density(1,i,:),sqrt(length(density)),sqrt(length(density(1,i,:))));
%     minDensity = (thisDensity>minDensityThreshold);    
%     
%     fig((i-1)*4 + 1) = figure((i-1)*4 + 1);hold off;
%     [~,h] = contourf(xq,yq,minDensity.*squeeze(fitq(1,i,:,:)),nlevels);
%     set(h,'LineColor','none');
%     c = colorbar;c.Label.String = 'Fitness';
%     caxis([0.8 1]);axis equal;axis([min(xlin) max(xlin) min(ylin) max(ylin)]);set(gca,'XTick',[]); set(gca,'YTick',[]);
%     ylabel(['Acquisition: ' strrep(xpNames{i},':GP','')]);
%     title('Rel. Bin Fitness');
%     
%     fig((i-1)*4 + 2) = figure((i-1)*4 + 2);hold off;
%     [~,h] = contourf(xq,yq,minDensity.*thisDensity,nlevels);
%     set(h,'LineColor','none')
%     c = colorbar;c.Label.String = 'pdf';
%     caxis([0 1]);axis equal;axis([min(xlin) max(xlin) min(ylin) max(ylin)]);set(gca,'XTick',[]); set(gca,'YTick',[]);
%     title('Acquired Samples');
%     
%     % Fitness and density predicted elites
%     fig((i-1)*4 + 3) = figure((i-1)*4 + 3);hold off;
%     thisDensity = reshape(density(2,i,:),sqrt(length(density)),sqrt(length(density(2,i,:))));
%     minDensity = (thisDensity>minDensityThreshold);    
%     [~,h] = contourf(xq,yq,minDensity.*thisDensity,nlevels);
%     set(h,'LineColor','none')
%     c = colorbar;c.Label.String = 'pdf';
%     caxis([0 1]);axis equal;axis([min(xlin) max(xlin) min(ylin) max(ylin)]);set(gca,'XTick',[]); set(gca,'YTick',[]);
%     title('Predicted Elites');
%     
%     
%     fig((i-1)*4 + 4) = figure((i-1)*4 + 4);hold off;
%     [~,h] = contourf(xq,yq,minDensity.*squeeze(fitq(2,i,:,:)),nlevels);
%     set(h,'LineColor','none');
%     c = colorbar;c.Label.String = 'Fitness';
%     caxis([0.8 1]);axis equal;axis([min(xlin) max(xlin) min(ylin) max(ylin)]);set(gca,'XTick',[]); set(gca,'YTick',[]);
%     ylabel(['Acquisition: ' strrep(xpNames{i},':GP','')]);
%     title('Rel. Bin Fitness');
%     
%     %fig((i-1)*3 + 3) = figure((i-1)*3 + 3);hold off;
%     %imagesc(cat(3,flipud(squeeze(binU(i,:,:))),flipud(squeeze(binV(i,:,:))),flipud(squeeze(binW(i,:,:)))));
%     %set(h,'LineColor','none');
%     %axis equal;axis tight;
%     %title('Bin IDs');
%     drawnow;
% end