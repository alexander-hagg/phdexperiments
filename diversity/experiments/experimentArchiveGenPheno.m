% EXPERIMENTGENETICNEUTRALITY
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2020; Last revision: 03-Feb-2020
%
%------------- BEGIN CODE --------------

clear;clc;
%%
DOF = 16;DOMAIN = 'npoly_ffd_imagePhenotype';
ALGORITHM = 'voronoi';
addpath(genpath('/home/alex/diversity'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);

d = domain(DOF);
p = defaultParamSet;
d.fitfun = d.fitfunPointSymm; % Multimodal function

p.nGens = 1024;

radialBoundAdaptation = 1;
axialBoundAdaptation = -1;
% Set domain ranges
disp(['Parameter Bounds: ' num2str(radialBoundAdaptation) ' / ' num2str(axialBoundAdaptation)]);
d.ranges(:,1) = [axialBoundAdaptation*ones(d.dof/2,1);-radialBoundAdaptation*pi*ones(d.dof/2,1)];
d.ranges(:,2) = [ ones(d.dof/2,1); radialBoundAdaptation*pi*ones(d.dof/2,1)];

maxBins = [10 25 50 100 200 400];
%maxBins = [10 100];
p.numInitSamples = maxBins;
nChildren = maxBins;

replicates = 5;

fname = ['archive_2_' int2str(p.maxBins(end)) '_DOF_' int2str(DOF) '.mat'];


%% Experiment
for rep=1:replicates
    disp(['Replicate ' int2str(rep) '/' int2str(replicates)]);
    % Initialize initial solution set for MAP-Elites
    sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = (rep-1)*p.numInitSamples + 1;
    initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
    [fitness,polygons] = d.fitfun(initSamples);
    for r = 1:length(maxBins)
        p.maxBins = maxBins(r);
        disp(['Replicate ' int2str(rep) '/' int2str(replicates)]);
        disp(['Max bins: ' int2str(p.maxBins)]);
        p.nChildren = p.maxBins;
        
        %% II) Illumination with QD (default)
        p.categorize = @(geno,pheno,p,d) categorize(pheno,d);
        features = p.categorize(initSamples, polygons, p, d);
        map = createMap(d, p,p.extraMapValues);
        [replaced, replacement] = nicheCompete(initSamples,fitness,map,d,p,features);
        map = updateMap(replaced,replacement,map,fitness,initSamples,features);
        [QDmap{r,rep}, percImproved, percValid, allMaps{r,rep}] = illuminate(map,d.fitfun,p,d);
        QD_genomes{r,rep} = reshape(QDmap{r,rep}.genes,[],d.dof); QD_genomes{r,rep}(all(isnan(QD_genomes{r,rep})'),:) = [];
        
        %%
        p.categorize = @(geno,pheno,p,d) geno;
        features = initSamples;
        map = createMap(d, p,p.extraMapValues);
        [replaced, replacement] = nicheCompete(initSamples,fitness,map,d,p,features);
        map = updateMap(replaced,replacement,map,fitness,initSamples,features);
        [QDGenomap{r,rep}, percImprovedGeno, percValidGeno, allMapsGeno{r,rep}] = illuminate(map,d.fitfun,p,d);
        
        QDGeno_genomes{r,rep} = reshape(QDGenomap{r,rep}.genes,[],d.dof); QDGeno_genomes{r,rep}(all(isnan(QDGeno_genomes{r,rep})'),:) = [];
        
    end
end

disp('Optimization Done');

%% Get auxiliary data
resolution = 64;
for rep=1:replicates
    disp(['Rep: ' int2str(rep)]);
    for r = 1:length(maxBins)
        disp(['R: ' int2str(r)]);
        [QD_fitness{r,rep},QD_polygons{r,rep},~] = d.fitfun(QD_genomes{r,rep});
        QD_features{r,rep} = categorize(QD_polygons{r,rep}, d)';
        [QD_phenotypes{r,rep}] = getPhenotypeBoolean(QD_polygons{r,rep},resolution);
        
        [QDGeno_fitness{r,rep},QDGeno_polygons{r,rep},~] = d.fitfun(QDGeno_genomes{r,rep});
        QDGeno_features{r,rep} = categorize(QDGeno_polygons{r,rep}, d)';
        [QDGeno_phenotypes{r,rep}] = getPhenotypeBoolean(QDGeno_polygons{r,rep},resolution);
        
    end
end

% Metrics
% Mean Distance to Nearest Neighbour
tic
for i=1:length(maxBins)
    for rep=1:replicates
        SOD_t(:,rep) = [metricSumOfDistancesNN(QDGeno_genomes{i,rep}); metricSumOfDistancesNN(QD_genomes{i,rep})];
        disp(['MODNN: ' num2str(toc)]);
    end
    if replicates==1
        MODNN(:,i) = (SOD_t');
    else
        MODNN(:,i) = median(SOD_t');
    end
end

tic
for i=1:length(maxBins)
    for rep=1:replicates
        SOD_t(:,rep) = [metricSumOfDistancesNN(QDGeno_phenotypes{i,rep},'hamming'); metricSumOfDistancesNN(QD_phenotypes{i,rep},'hamming')];
        disp(['MODNNp: ' num2str(toc)]);
    end
    if replicates==1
        MODNNp(:,i) = (SOD_t');
    else
        MODNNp(:,i) = median(SOD_t');
    end
end


% Solow-Polasky Diversity
clear SPD SPDp;
tic
theta = 1; %1 0.5
for i=1:length(maxBins)
    for rep=1:replicates
        SPD_t(:,rep) = [metricSPD(QDGeno_genomes{i,rep},'euclidean',theta); metricSPD(QD_genomes{i,rep},'euclidean',theta)];
        disp(['SPD: ' num2str(toc)]);
    end
    if replicates==1
        SPD(:,i) = (SPD_t');
    else
        SPD(:,i) = median(SPD_t');
    end
end

tic
theta = 100; %100 50
for i=1:length(maxBins)
    for rep=1:replicates
        SPD_t(:,rep) = [metricSPD(QDGeno_phenotypes{i,rep},'hamming',theta); metricSPD(QD_phenotypes{i,rep},'hamming',theta)];
        disp(['SPDp: ' num2str(toc)]);
    end
    if replicates==1
        SPDp(:,i) = (SPD_t');
    else
        SPDp(:,i) = median(SPD_t');
    end
end

% Pure Diversity
for i=1:length(maxBins)
    for rep=1:replicates
        PD_t(:,rep) = [metricPD(QDGeno_genomes{i,rep}); metricPD(QD_genomes{i,rep})];
        disp(['PD: ' num2str(toc)]);
    end
    if replicates==1
        PD(:,i) = (PD_t');
    else
        PD(:,i) = median(PD_t');
    end
end

for i=1:length(maxBins)
    for rep=1:replicates
        PD_t(:,rep) = [metricPD(QDGeno_phenotypes{i,rep},'hamming'); metricPD(QD_phenotypes{i,rep},'hamming')];
        disp(['PDp: ' num2str(toc)]);
    end
    if replicates==1
        PDp(:,i) = (PD_t');
    else
        PDp(:,i) = median(PD_t');
    end
end

%% Fitness
clear fit fit_t;
for i=1:length(maxBins)
    fit(:,i) =[median([QDGeno_fitness{i,:}]); median([QD_fitness{i,:}])];
end

%%
save(fname);
system(['scp -r ' fname ' ahagg2s@wr0.wr.inf.h-brs.de:/home/ahagg2s/diversity']);

%% Plots
set(0,'DefaultFigureWindowStyle','default')
clear figs;
lw = 2;

% Mean Distance NN
figs(1) = figure(1);hold off;
semilogx(repmat(maxBins(2:end),2,1)',MODNN(:,2:end)', 'LineWidth', lw);
grid on; %legend('MAP-Elites Geno','MAP-Elites Pheno', 'Location','NorthWest');
title('Sum Genetic Distance to Nearest Neighbour');
xlabel('Max. # Bins');
axis([maxBins(2) maxBins(end) 0 2000]);
ax = gca; ax.XTick = maxBins(2:end);

figs(end+1) = figure;hold off;
semilogx(repmat(maxBins(2:end),2,1)',MODNNp(:,2:end)', 'LineWidth', lw);
grid on;% legend('MAP-Elites Geno','MAP-Elites Pheno', 'Location','NorthWest');
title('Sum Phenotypic Distance to Nearest Neighbour');
xlabel('Max. # Bins');
axis([maxBins(2) maxBins(end) 0 15]);
ax = gca; ax.XTick = maxBins(2:end);

% Solow-Polasky Diversity
figs(end+1) = figure;hold off;
semilogx(repmat(maxBins(2:end),2,1)',SPD(:,2:end)', 'LineWidth', lw);
grid on;% legend('MAP-Elites Geno','MAP-Elites Pheno', 'Location','NorthWest');
title('Genetic Solow-Polasky Diversity');
xlabel('Max. # Bins');
axis([maxBins(2) maxBins(end) 0 350]);
ax = gca; ax.XTick = maxBins(2:end);

figs(end+1) = figure;hold off;
semilogx(repmat(maxBins(2:end),2,1)',SPDp(:,2:end)', 'LineWidth', lw);
grid on; %legend('MAP-Elites Geno','MAP-Elites Pheno', 'Location','NorthWest');
title('Phenotypic Solow-Polasky Diversity');
xlabel('Max. # Bins');
axis([maxBins(2) maxBins(end) 0 300]);
ax = gca; ax.XTick = maxBins(2:end);

% Pure Diversity
figs(end+1) = figure;hold off;
semilogx(repmat(maxBins(2:end),2,1)',PD(:,2:end)', 'LineWidth', lw);
grid on; %legend('MAP-Elites Geno','MAP-Elites Pheno', 'Location','NorthWest');
title('Genetic Pure Diversity');
xlabel('Max. # Bins');
axis([maxBins(2) maxBins(end) 0 15e13]);
ax = gca; ax.XTick = maxBins(2:end);

figs(end+1) = figure;hold off;
semilogx(repmat(maxBins(2:end),2,1)',PDp(:,2:end)', 'LineWidth', lw);
grid on; %legend('MAP-Elites Geno','MAP-Elites Pheno', 'Location','NorthWest');
title('Phenotypic Pure Diversity');
xlabel('Max. # Bins');
axis([maxBins(2) maxBins(end) 0 15]);
ax = gca; ax.XTick = maxBins(2:end);

% Fitness
figs(end+1) = figure;hold off;
semilogx(repmat(maxBins(2:end),2,1)',fit(:,2:end)', 'LineWidth', lw);
grid on; %legend('MAP-Elites Geno','MAP-Elites Pheno', 'Location','NorthWest');
title('Avg. Fitness');
xlabel('Max. # Bins');
axis([maxBins(2) maxBins(end) 0.9 1.05]);
ax = gca; ax.XTick = maxBins(2:end);

save_figures(figs,'.','QDgenopheno',11,[4 3])

%% Fitness
rep = 2
for r=1:length(maxBins)
    figs(end+1) = figure;hold off;
    plot(QDGenomap{r,rep}.stats.fitnessMean);
    hold on;
    plot(QDmap{r,rep}.stats.fitnessMean);
    title(['Mean Fitness Gen']);
    xlabel('Gens'); ylabel('Mean Fitness (%)');
    legend('Genetic Map','Phenotypic Map','Location','SouthEast');
    ax = gca;
    ax.YLim = [0 1];
    grid on;
end

%% Show examples
rep=1;
for i=1:length(maxBins)
    
    figs(end+1) = figure;hold off;
    showPhenotype(QDGeno_genomes{i,rep},d,figs(end));
    title(['VE Genetic, # solutions: ' num2str(maxBins(i))]);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
    figs(end+1) = figure;hold off;
    showPhenotype(QD_genomes{i,rep},d,figs(end));
    title(['VE Phenotypical, # solutions: ' num2str(maxBins(i))]);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
end

save_figures(figs,'.','QDgenopheno',11,[4 3])

%% Maps
set(0,'DefaultFigureWindowStyle','default')
jj = 1;
for ii=1:length(maxBins)
    
    figs((ii-1)*2+1) = figure((ii-1)*2+1);ax = gca;
    [figHandle, imageHandle, cHandle] = viewMap(QDmap{r,rep},'fitness',d,ax);caxis([0.5 1]);
    cHandle.Label.String = ['Symmetry\rightarrow'];
    title('Phenotype Diversity');
    
    figs((ii-1)*2+2) = figure((ii-1)*2+2);ax = gca;
    [figHandle, imageHandle, cHandle] = viewMap(QDGenomap{r,rep},'fitness',d,ax);caxis([0.5 1]);
    cHandle.Label.String = ['Symmetry\rightarrow'];
    title('Genotype Diversity');
    
end

%%
save_figures(figs,'.','QDgenopheno_maps',11,[6 5])



%% Genetic Space Reach
clear figs;
cmap = [0 0 0; 1.0 0.2 0.2; 0 0 0];
lineStyle = {':','-','--'};
figs(1) = figure(1);
for i=1:size(QDGeno_genomes,1)
    clear mins maxs
    for rep=1:replicates
        mins(rep,:,:) = [min(QDGeno_genomes{i,rep});min(QD_genomes{i,rep})];
        maxs(rep,:,:) = [max(QDGeno_genomes{i,rep});max(QD_genomes{i,rep})];
    end
    
    mins = squeeze(median(mins,1));
    maxs = squeeze(median(maxs,1));
    
    subplot(size(QDGeno_genomes,1),1,i);hold off;
    for j=1:2
        h(j) = plot(mins(j,:),'LineStyle',lineStyle{j},'Color',cmap(j,:));
        hold on;
        plot(maxs(j,:),'LineStyle',lineStyle{j},'Color',cmap(j,:));
    end
    grid on;
    ax = gca;
    ax.XLim = [0 20];
    ax.XTick = 1:DOF;
    
    ticks = [axialBoundAdaptation, 1, -radialBoundAdaptation*pi,radialBoundAdaptation*pi];
    ax.YLim = [min(ticks) - 0.5, max(ticks) + 0.5];
    %ax.YTick = unique(sort(ticks,'ascend'));
    legend(h, 'QD geno','QD pheno');
    if i==1
        title('Genetic Parameters Covered');
    end
    
end

save_figures(figs,'.','QDgenopheno_parameters',11,[6 5])
