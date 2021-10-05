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

%%

d = domain(DOF);
p = defaultParamSet;
d.fitfun = d.fitfunPointSymm; % Multimodal function

p.nGens = 1024;
p.maxBins = 400;
p.numInitSamples = p.maxBins;
p.nChildren = p.numInitSamples;

axialBoundAdaptations =   [0.0   0.0   -0.25  -0.50 -1.00];
radialBoundAdaptations =  [0.05  0.125  0.25   0.50  1.00];
p.categorize = @(g,p,d) categorize(p,d);
fname = ['neutrality' int2str(p.maxBins) '_DOF_' int2str(DOF) '.mat'];

replicates = 5;

%%
for rep=1:replicates
    for r = 1:length(radialBoundAdaptations)
        disp(['Replicate ' int2str(rep) '/' int2str(replicates)]);
        % Set domain ranges
        radialBoundAdaptation = radialBoundAdaptations(r);
        axialBoundAdaptation = axialBoundAdaptations(r);
        disp(['Parameter Bounds: ' num2str(radialBoundAdaptation) ' / ' num2str(axialBoundAdaptation)]);
        d.ranges(:,1) = [ axialBoundAdaptation*ones(d.dof/2,1);-radialBoundAdaptation*pi*ones(d.dof/2,1)];
        d.ranges(:,2) = [ ones(d.dof/2,1); radialBoundAdaptation*pi*ones(d.dof/2,1)];
        disp(['Feature bounds: ' mat2str(d.featureMin) ' - ' mat2str(d.featureMax)]);
        
        
        % Reinitialize initial solution set for NSGA2 and VE
        sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = (rep-1)*p.numInitSamples + 1;
        initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
        
        [fitness,polygons] = d.fitfun(initSamples);
        features = p.categorize(initSamples,polygons, d);
        
        %% I) Multiobjective optimization with NSGA-II
        numObjectives = evaluate_objective;
        pc = 0.9; pm = 1/d.dof;
        nsgaFun = @(x) evaluate_objective(x,d,false);
        [R, Rfit, Rrank, stats] = NSGAII(p.nChildren,p.nGens,pc,pm,p.mutSigma,nsgaFun,d.dof,d.ranges(:,1),d.ranges(:,2),initSamples);
        NSGA_genomes{r,rep} = R;
        
        %% II) Illumination with QD (default)
        p.categorize = @(g,p,d) categorize(p,d);
        map = createMap(d, p,p.extraMapValues);
        [replaced, replacement] = nicheCompete(initSamples,fitness,map,d,p,features);
        map = updateMap(replaced,replacement,map,fitness,initSamples,features);
        [QDmap, percImproved, percValid, allMaps, ~, statsQD{r,rep}] = illuminate(map,d.fitfun,p,d);
        QD_genomes{r,rep} = reshape(QDmap.genes,[],d.dof); QD_genomes{r,rep}(all(isnan(QD_genomes{r,rep})'),:) = [];
        
        %% III) Niching in Genetic Space (Python)
        numOuterIters = p.maxBins;
        numInnerEvals = ceil(p.nChildren*p.nGens/numOuterIters);
        localSearchMethod = 'L-BFGS-B'; %'L-BFGS-B'; 'COBYLA'
        delete '~/diversity/MMO/basinhopping/*.csv';
        cd('~/diversity/MMO/basinhopping');
        system(['python3 RLS.py -n ' int2str(numInnerEvals) ' -i ' int2str(numOuterIters) ' -l ' localSearchMethod ' -a ' num2str(axialBoundAdaptation) ' -m ' num2str(radialBoundAdaptation) ' -d ' int2str(DOF)]);
        cd('~/diversity');
        RLS_genomes{r,rep} = csvread('MMO/basinhopping/optima.csv');
        %
    end
end

disp('Optimization Done');

%% Get auxiliary data
resolution = 64;
for rep=1:replicates
    disp(['Rep: ' int2str(rep)]);
    for r = 1:length(radialBoundAdaptations)
        disp(['R: ' int2str(r)]);
        [NSGA_fitness{r,rep},NSGA_polygons{r,rep},~] = d.fitfun(NSGA_genomes{r,rep});
        NSGA_features{r,rep} = categorize(NSGA_polygons{r,rep}, d)';
        [NSGA_phenotypes{r,rep}] = getPhenotypeBoolean(NSGA_polygons{r,rep},resolution);
        
        [QD_fitness{r,rep},QD_polygons{r,rep},~] = d.fitfun(QD_genomes{r,rep});
        QD_features{r,rep} = categorize(QD_polygons{r,rep}, d)';
        [QD_phenotypes{r,rep}] = getPhenotypeBoolean(QD_polygons{r,rep},resolution);
        
        [RLS_fitness{r,rep},RLS_polygons{r,rep}] = d.fitfun(RLS_genomes{r,rep});
        RLS_features{r,rep} = categorize(RLS_polygons{r,rep}, d);
        [RLS_phenotypes{r,rep}] = getPhenotypeBoolean(RLS_polygons{r,rep},resolution);
    end
end

%% Metrics
% Mean Distance to Nearest Neighbour
clear MODNN MODNNp SPD SPDp PD PDp SOD_t
tic
for i=1:length(radialBoundAdaptations)
    for rep=1:replicates
        SOD_t(:,rep) = [metricSumOfDistancesNN(NSGA_genomes{i,rep}); metricSumOfDistancesNN(QD_genomes{i,rep}); metricSumOfDistancesNN(RLS_genomes{i,rep})];
        disp(['MODNN: ' num2str(toc)]);
    end
    if replicates==1
        MODNN(:,i) = (SOD_t');
    else
        MODNN(:,i) = median(SOD_t');
    end
end

tic
for i=1:length(radialBoundAdaptations)
    for rep=1:replicates
        SOD_t(:,rep) = [metricSumOfDistancesNN(NSGA_phenotypes{i,rep},'hamming'); metricSumOfDistancesNN(QD_phenotypes{i,rep},'hamming'); metricSumOfDistancesNN(RLS_phenotypes{i,rep},'hamming')];
        disp(['MODNNp: ' num2str(toc)]);
    end
    if replicates==1
        MODNNp(:,i) = (SOD_t');
    else
        MODNNp(:,i) = median(SOD_t');
    end
end

%% Solow-Polasky Diversity
theta = 1;
tic
for i=1:length(radialBoundAdaptations)
    for rep=1:replicates
        SPD_t(:,rep) = [metricSPD(NSGA_genomes{i,rep},'euclidean',theta); metricSPD(QD_genomes{i,rep},'euclidean',theta); metricSPD(RLS_genomes{i,rep},'euclidean',theta)];
        disp(['SPD: ' num2str(toc)]);
    end
    if replicates==1
        SPD(:,i) = (SPD_t');
    else
        SPD(:,i) = median(SPD_t');
    end
end

theta = 100;
tic
for i=1:length(radialBoundAdaptations)
    for rep=1:replicates
        SPD_t(:,rep) = [metricSPD(NSGA_phenotypes{i,rep},'hamming',theta); metricSPD(QD_phenotypes{i,rep},'hamming',theta); metricSPD(RLS_phenotypes{i,rep},'hamming',theta)];
        disp(['SPDp: ' num2str(toc)]);
    end
    if replicates==1
        SPDp(:,i) = (SPD_t');
    else
        SPDp(:,i) = median(SPD_t');
    end
end

% Pure Diversity

tic;
for i=1:length(radialBoundAdaptations)
    for rep=1:replicates
        SPD_t(:,rep) = [metricPD(NSGA_genomes{i,rep}); metricPD(QD_genomes{i,rep}); metricPD(RLS_genomes{i,rep})];
        disp(['PD: ' num2str(toc)]);
    end
    if replicates==1
        PD(:,i) = (SPD_t');
    else
        PD(:,i) = median(SPD_t');
    end
end

tic
for i=1:length(radialBoundAdaptations)
    for rep=1:replicates
        SPD_t(:,rep) = [metricPD(NSGA_phenotypes{i,rep},'hamming'); metricPD(QD_phenotypes{i,rep},'hamming'); metricPD(RLS_phenotypes{i,rep},'hamming')];
        disp(['PDp: ' num2str(toc)]);
    end
    if replicates==1
        PDp(:,i) = (SPD_t');
    else
        PDp(:,i) = median(SPD_t');
    end
end



%% Fitness
clear fitSum fitMean fitMedian fitStd;
for r=1:length(radialBoundAdaptations)
    fitMedian(r,:) = [median([NSGA_fitness{r,:}]),median([QD_fitness{r,:}]),median([RLS_fitness{r,:}])];
end
%fitSum = squeeze(mean(fitSum,2));
%fitMean = squeeze(mean(fitMean,2));
fitMedian = fitMedian;
%fitStd = squeeze(mean(fitStd,2));


%%
save(fname);
system(['scp -r ' fname ' ahagg2s@wr0.wr.inf.h-brs.de:/home/ahagg2s/diversity']);

%% Plots
set(0,'DefaultFigureWindowStyle','default')
clear figs;
lw = 2;

% Sum Distance NN
figs(3) = figure(3);hold off;
semilogx(repmat(radialBoundAdaptations,3,1)',MODNN', 'LineWidth', lw);
hold on;
plot([0.125,0.125],[0 2500],'k--', 'LineWidth', lw);
grid on; %legend('NSGA2','VE',['RLS-' localSearchMethod],'Neutrality','Location','NorthWest'); 
title('Sum Genetic Distance to Nearest Neighbour');
xlabel('Max. Radial Deviation (Genetic Neutrality)');% ylabel('Euclidean Distance');
axis([0 radialBoundAdaptations(end) 0 2000]);
ax = gca; ax.XTick = radialBoundAdaptations;

figs(4) = figure(4);hold off;
semilogx(repmat(radialBoundAdaptations,3,1)',MODNNp', 'LineWidth', lw);
hold on;
plot([0.125,0.125],[0 2500],'k--', 'LineWidth', lw);
grid on; %legend('NSGA2','VE',['RLS-' localSearchMethod],'Neutrality','Location','NorthWest');
title('Sum Phenotypic Distance to Nearest Neighbour');
xlabel('Max. Radial Deviation (Genetic Neutrality)');%ylabel('Hamming Distance');
axis([0 radialBoundAdaptations(end) 0 15]);
ax = gca; ax.XTick = radialBoundAdaptations;

% Solow-Polasky Diversity
figs(5) = figure(5);hold off;
semilogx(repmat(radialBoundAdaptations,3,1)',SPD', 'LineWidth', lw);
hold on;
plot([0.125,0.125],[0 2500],'k--', 'LineWidth', lw);
grid on; %legend('NSGA2','VE',['RLS-' localSearchMethod],'Neutrality','Location','NorthWest');
title('Genetic Solow-Polasky Diversity');
xlabel('Max. Radial Deviation (Genetic Neutrality)'); %ylabel('Euclidean Distance');
axis([0 radialBoundAdaptations(end) 0 350]);
ax = gca; ax.XTick = radialBoundAdaptations;

figs(6) = figure(6);hold off;
semilogx(repmat(radialBoundAdaptations,3,1)',SPDp', 'LineWidth', lw);
hold on;
plot([0.125,0.125],[0 2500],'k--', 'LineWidth', lw);
grid on;% legend('NSGA2','VE',['RLS-' localSearchMethod],'Neutrality','Location','NorthWest');
title('Phenotypic Solow-Polasky Diversity');
xlabel('Max. Radial Deviation (Genetic Neutrality)');%ylabel('Hamming Distance');
axis([0 radialBoundAdaptations(end) 0 300]);
ax = gca; ax.XTick = radialBoundAdaptations;

% Pure Diversity
figs(7) = figure(7);hold off;
semilogx(repmat(radialBoundAdaptations,3,1)',PD', 'LineWidth', lw);
hold on;
plot([0.125,0.125],[0 3e14],'k--', 'LineWidth', lw);
grid on;% legend('NSGA2','VE',['RLS-' localSearchMethod],'Neutrality','Location','NorthWest');
title('Genetic Pure Diversity');
xlabel('Max. Radial Deviation (Genetic Neutrality)');% ylabel('Minkowski Distance');
axis([0 radialBoundAdaptations(end) 0 15e13]);
ax = gca; ax.XTick = radialBoundAdaptations;

figs(8) = figure(8);hold off;
semilogx(repmat(radialBoundAdaptations,3,1)',PDp', 'LineWidth', lw);
hold on;
plot([0.125,0.125],[0 2500],'k--', 'LineWidth', lw);
grid on; %legend('NSGA2','VE',['RLS-' localSearchMethod],'Neutrality','Location','NorthWest'); 
title('Phenotypic Pure Diversity');
xlabel('Max. Radial Deviation (Genetic Neutrality)');%ylabel('Hamming Distance');
axis([0 radialBoundAdaptations(end) 0 15]);
ax = gca; ax.XTick = radialBoundAdaptations;

% Fitness
figs(9) = figure(9);hold off;
semilogx(repmat(radialBoundAdaptations,3,1)',fitMedian, 'LineWidth', lw);
hold on;
plot([0.125,0.125],[0 2500],'k--', 'LineWidth', lw);
title(['median Fitness']);
xlabel('Max. Radial Deviation (Genetic Neutrality)');% ylabel('Fitness');
grid on; %legend('NSGA2','VE',['RLS-' localSearchMethod],'Neutrality','Location','SouthWest');
axis([0 radialBoundAdaptations(end) 0.9 1.05]);
ax = gca; ax.XTick = radialBoundAdaptations;

% figs(10) = figure(10);hold off;
% semilogx(repmat(radialBoundAdaptations,3,1)',fitStd, 'LineWidth', lw);
% hold on;
% plot([0.125,0.125],[0 2500],'k--', 'LineWidth', lw);
% title(['\sigma Fitness']);
% xlabel('Max. Radial Deviation (Genetic Neutrality)'); %ylabel('Fitness');
% grid on;% legend('NSGA2','VE',['RLS-' localSearchMethod],'Neutrality','Location','NorthWest');
% axis([0 radialBoundAdaptations(end) 0 0.05]);
% ax = gca; ax.XTick = radialBoundAdaptations;

%% Show examples
rep=1;
showing = [2 3 5];
for bla=1:length(showing)
    i=showing(bla);
    %for i=1:length(radialBoundAdaptations)
    
    %varThreshold = 95; [~,score] = metricCompressibilityPCA([NSGA_phenotypes{i,rep};NSGA_phenotypes_NOX{i,rep};QD_phenotypes{i,rep};RLS_phenotypes{i,rep}],varThreshold);score = score*2;
    %score = tsne([NSGA_phenotypes{i,rep};NSGA_phenotypes_NOX{i,rep};QD_phenotypes{i,rep};RLS_phenotypes{i,rep}],'Algorithm','barneshut','Standardize',true,'Perplexity',30,'Theta',0.2);score = score*0.5;
    
    %axRanges = [min(score(:))*1.05 max(score(:))*1.05 min(score(:))*1.05 max(score(:))*1.05];
    
    szNSGA = size(NSGA_phenotypes{i,rep},1);
    szQD = size(QD_phenotypes{i,rep},1);
    szRLS = size(RLS_phenotypes{i,rep},1);
    sel1 = 1;
    sel2 = 2;
    
    figs(end+1) = figure;hold off;
    %placement = score(1:szNSGA,[sel1,sel2]);
    showPhenotype(NSGA_genomes{i,rep},d,figs(end));
    %axis(axRanges);
    title(['NSGA2, radial multiplier: ' num2str(radialBoundAdaptations(i))]);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
    figs(end+1) = figure;hold off;
    %placement = score(szNSGA+szNSGAinv+1:szNSGA+szNSGAinv+szQD,[sel1,sel2]);
    showPhenotype(QD_genomes{i,rep},d,figs(end));
    title(['VE, radial multiplier: ' num2str(radialBoundAdaptations(i))]);
    %axis(axRanges);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
    figs(end+1) = figure;hold off;
    %placement = score(szNSGA+szNSGAinv+szQD+1:end,[sel1,sel2]);
    showPhenotype(RLS_genomes{i,rep},d,figs(end));
    title([['RLS-' localSearchMethod] ', radial multiplier: ' num2str(radialBoundAdaptations(i))]);
    %axis(axRanges);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
end

%%

save_figures(figs,'.','3algs',11,[4 3])


%% Genetic Space Reach
clear figs;
cmap = [0 0 0; 1.0 0.2 0.2; 0 0 0];
lineStyle = {':','-','--'};
figs(1) = figure(1);
for i=1:size(NSGA_genomes,1)
    clear mins maxs
    for rep=1:replicates
        mins(rep,:,:) = [min(NSGA_genomes{i,rep});min(QD_genomes{i,rep});min(RLS_genomes{i,rep})];
        maxs(rep,:,:) = [max(NSGA_genomes{i,rep});max(QD_genomes{i,rep});max(RLS_genomes{i,rep})];
    end
    
    mins = squeeze(median(mins,1));
    maxs = squeeze(median(maxs,1));
    
    subplot(size(NSGA_genomes,1),1,i);hold off;
    for j=1:3
        h(j) = plot(mins(j,:),'LineStyle',lineStyle{j},'Color',cmap(j,:));
        hold on;
        plot(maxs(j,:),'LineStyle',lineStyle{j},'Color',cmap(j,:));
    end
    grid on;
    ax = gca;
    ax.XLim = [0 20];
    ax.XTick = 1:DOF;
    
    ticks = [axialBoundAdaptations(i), 1, -radialBoundAdaptations(i)*pi,radialBoundAdaptations(i)*pi];
    ax.YLim = [min(ticks) - 0.5, max(ticks) + 0.5];
    legend(h, 'NSGA-II','QD','RLS');
    if i==1
        title('Genetic Parameters Covered');
    end
    
end

save_figures(figs,'.','3algs_reachGeneticSpace',11,[8 8])



%% Maps

set(0,'DefaultFigureWindowStyle','default')

clear figs
for ii = 1:5
    jj = 1;
    
    NSGAmap = createMap(d, p);
    [replaced, replacement, ~]              = nicheCompete(NSGA_genomes{ii,jj}, NSGA_fitness{ii,jj}, NSGAmap, d, p, NSGA_features{ii,jj}');
    NSGAmap                                 = updateMap(replaced,replacement,NSGAmap,NSGA_fitness{ii,jj}',NSGA_genomes{ii,jj},NSGA_features{ii,jj}');
    
    QDmap = createMap(d, p);
    [replaced, replacement, ~] = nicheCompete(QD_genomes{ii,jj}, QD_fitness{ii,jj}, QDmap, d, p, QD_features{ii,jj}');
    QDmap                      = updateMap(replaced,replacement,QDmap,QD_fitness{ii,jj}',QD_genomes{ii,jj},QD_features{ii,jj}');
    
    RLSmap = createMap(d, p);
    [replaced, replacement, ~] = nicheCompete(RLS_genomes{ii,jj}, RLS_fitness{ii,jj}, RLSmap, d, p, RLS_features{ii,jj});
    RLSmap                     = updateMap(replaced,replacement,RLSmap,RLS_fitness{ii,jj}',RLS_genomes{ii,jj},RLS_features{ii,jj});
    
    [figHandle, imageHandle, cHandle] = viewMap(NSGAmap,'fitness', d);caxis([0.5 1]);
    cHandle.Label.String = ['Symmetry\rightarrow'];
    figs((ii-1)*3+1) = gcf;
    
    [figHandle, imageHandle, cHandle] = viewMap(QDmap,'fitness', d);caxis([0.5 1]);
    cHandle.Label.String = ['Symmetry\rightarrow'];
    figs((ii-1)*3+2) = gcf;
    
    [figHandle, imageHandle, cHandle] = viewMap(RLSmap,'fitness', d);caxis([0.5 1]);
    cHandle.Label.String = ['Symmetry\rightarrow'];
    figs((ii-1)*3+3) = gcf;
    
end

save_figures(figs,'.','3algs_maps',12,[4 4])




