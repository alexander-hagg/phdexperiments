% EXPERIMENTAUTOARCHIVE
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Mar 2020; Last revision: 13-Mar-2020
%
%------------- BEGIN CODE --------------

clear;clc;
%%
DOF = 16;DOMAIN = 'npoly_ffd_imagePhenotype';
ALGORITHM = 'voronoi';
addpath(genpath('/home/alex/diversity'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);

d = domain(DOF);
pStd = defaultParamSet;

d.fitfun = d.fitfunPointSymm; % Multimodal function

radialBoundAdaptation = 1; axialBoundAdaptation = -1;
% Set domain ranges
disp(['Parameter Bounds: ' num2str(radialBoundAdaptation) ' / ' num2str(axialBoundAdaptation)]);
d.ranges(:,1) = [axialBoundAdaptation*ones(d.dof/2,1);-radialBoundAdaptation*pi*ones(d.dof/2,1)];
d.ranges(:,2) = [ ones(d.dof/2,1); radialBoundAdaptation*pi*ones(d.dof/2,1)];

replicates = 5;
fname = ['autov5' int2str(pStd.maxBins(end)) '_DOF_' int2str(DOF) '.mat'];
maxBins = [25 50 100 200 400];

nGens = 512;
numIterations = 2;

%% Experiment
for r=1:length(maxBins)
    pStd.nGens = nGens;
    pStd.maxBins = maxBins(r); pStd.numInitSamples = pStd.maxBins; pStd.nChildren = pStd.maxBins;
    
    % Reset Manual QD nGen to poem nIter * nGen
    p = pStd;
    p.nGens = nGens * numIterations;
    
    %%
    for rep=1:replicates
        disp(['Replicate ' int2str(rep) '/' int2str(replicates)]);
        % Initialize initial solution set for MAP-Elites
        sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = (rep-1)*p.numInitSamples + 1;
        initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
        [fitness,polygons] = d.fitfun(initSamples);
        
        %% I) Illumination with QD (default)
        p.categorize = @(geno,pheno,p,d) categorize(pheno,d);
        features = p.categorize(initSamples,polygons, p, d);
        map = createMap(d, p,p.extraMapValues);
        [replaced, replacement] = nicheCompete(initSamples,fitness,map,d,p,features);
        map = updateMap(replaced,replacement,map,fitness,initSamples,features);
        [MANUmap{r,rep}, percImproved, percValid, allMaps{r,rep}] = illuminate(map,d.fitfun,p,d);
        MANU_genomes{r,rep} = MANUmap{r,rep}.genes;
        [MANU_fitness{r,rep},MANU_polygons{r,rep},~] = d.fitfun(MANU_genomes{r,rep});
        MANU_features{r,rep} = MANUmap{r,rep}.features;
        MANU_phenotypes{r,rep} = getPhenotypeBoolean(MANU_polygons{r,rep},d.resolution);
        
        %% II) Illumination with QD (poem)
        m = cfgLatentModel('data/workdir',d.resolution);
        poemCfg = poemParamSet(pStd,m);
        poemCfg.numInitSamples                = pStd.numInitSamples;
        poemCfg.numIterations                 = numIterations;
        poemCfg.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);
        
        [AUTOmap{r,rep},configs{r,rep},stats{r,rep}] = poem(initSamples,polygons,fitness,poemCfg,d);
        AUTO_genomes{r,rep} = AUTOmap{r,rep}.genes;
        [AUTO_fitness{r,rep},AUTO_polygons{r,rep},~] = d.fitfun(AUTO_genomes{r,rep});
        AUTO_features{r,rep} = AUTOmap{r,rep}.features;
        AUTO_phenotypes{r,rep} = getPhenotypeBoolean(AUTO_polygons{r,rep},d.resolution);
        
        %% lat dim 5
        m = getAEConfig('data/workdir',d.resolution,5);
        poemCfg = poemParamSet(pStd,m);
        poemCfg.numInitSamples                = pStd.numInitSamples;
        poemCfg.numIterations                 = numIterations;
        poemCfg.categorize = @(geno,pheno,p,d) predictAEfeatures(pheno,p.model);
        
        [AUTOmap5{r,rep},configs5{r,rep},stats5{r,rep}] = poem(initSamples,polygons,fitness,poemCfg,d);
        AUTO_genomes5{r,rep} = AUTOmap5{r,rep}.genes;
        [AUTO_fitness5{r,rep},AUTO_polygons5{r,rep},~] = d.fitfun(AUTO_genomes5{r,rep});
        AUTO_features5{r,rep} = AUTOmap5{r,rep}.features;
        AUTO_phenotypes5{r,rep} = getPhenotypeBoolean(AUTO_polygons5{r,rep},d.resolution);
        
        %% lat dim 10
        m = getAEConfig('data/workdir',d.resolution,10);
        poemCfg = poemParamSet(pStd,m);
        poemCfg.numInitSamples                = pStd.numInitSamples;
        poemCfg.numIterations                 = numIterations;
        poemCfg.categorize = @(geno,pheno,p,d) predictAEfeatures(pheno,p.model);
        
        [AUTOmap10{r,rep},configs10{r,rep},stats10{r,rep}] = poem(initSamples,polygons,fitness,poemCfg,d);
        AUTO_genomes10{r,rep} = AUTOmap10{r,rep}.genes;
        [AUTO_fitness10{r,rep},AUTO_polygons10{r,rep},~] = d.fitfun(AUTO_genomes10{r,rep});
        AUTO_features10{r,rep} = AUTOmap10{r,rep}.features;
        AUTO_phenotypes10{r,rep} = getPhenotypeBoolean(AUTO_polygons10{r,rep},d.resolution);
        
    end
end

save(fname);
%system(['scp -r ' fname ' ahagg2s@wr0.wr.inf.h-brs.de:/home/ahagg2s/diversity']);

disp('Optimization Done');

%% Analysis
clear MODNN MODNNp SPD SPDp PD PDp SPD_t
% Sum Distance to Nearest Neighbour
tic
for i=1:length(maxBins)
    for rep=1:replicates
        SOD_t(:,rep) = [metricSumOfDistancesNN(MANU_genomes{i,rep}); metricSumOfDistancesNN(AUTO_genomes{i,rep}); metricSumOfDistancesNN(AUTO_genomes5{i,rep}); metricSumOfDistancesNN(AUTO_genomes10{i,rep})];
        disp(['MODNN: ' num2str(toc)]);
    end
    if replicates==1
        MODNN(:,i) = (SOD_t');
    else
        MODNN(:,i) = median(SOD_t');
    end
end
clear SOD_t

tic
for i=1:length(maxBins)
    for rep=1:replicates
        SOD_t(:,rep) = [metricSumOfDistancesNN(MANU_phenotypes{i,rep},'hamming'); metricSumOfDistancesNN(AUTO_phenotypes{i,rep},'hamming'); metricSumOfDistancesNN(AUTO_phenotypes5{i,rep},'hamming'); metricSumOfDistancesNN(AUTO_phenotypes10{i,rep},'hamming')];
        disp(['MODNNp: ' num2str(toc)]);
    end
    if replicates==1
        MODNNp(:,i) = (SOD_t');
    else
        MODNNp(:,i) = median(SOD_t');
    end
end
clear SOD_t

% Solow-Polasky Diversity
theta = 1;
tic
for i=1:length(maxBins)
    for rep=1:replicates
        SPD_t(:,rep) = [metricSPD(MANU_genomes{i,rep},'euclidean',theta); metricSPD(AUTO_genomes{i,rep},'euclidean',theta); metricSPD(AUTO_genomes5{i,rep},'euclidean',theta); metricSPD(AUTO_genomes10{i,rep},'euclidean',theta)];
        disp(['SPD: ' num2str(toc)]);
    end
    if replicates==1
        SPD(:,i) = (SPD_t');
    else
        SPD(:,i) = median(SPD_t');
    end
end
clear SPD_t


theta = 100;
tic
for i=1:length(maxBins)
    for rep=1:replicates
        SPD_t(:,rep) = [metricSPD(MANU_phenotypes{i,rep},'hamming', theta); metricSPD(AUTO_phenotypes{i,rep},'hamming', theta); metricSPD(AUTO_phenotypes5{i,rep},'hamming', theta); metricSPD(AUTO_phenotypes10{i,rep},'hamming', theta)];
        disp(['SPDp: ' num2str(toc)]);
    end
    if replicates==1
        SPDp(:,i) = (SPD_t');
    else
        SPDp(:,i) = median(SPD_t');
    end
end

clear SPD_t;

% Pure Diversity

tic;
for i=1:length(maxBins)
    for rep=1:replicates
        SPD_t(:,rep) = [metricPD(MANU_genomes{i,rep}); metricPD(AUTO_genomes{i,rep}); metricPD(AUTO_genomes5{i,rep}); metricPD(AUTO_genomes10{i,rep})];
        disp(['PD: ' num2str(toc)]);
    end
    if replicates==1
        PD(:,i) = (SPD_t');
    else
        PD(:,i) = median(SPD_t');
    end
end

clear SPD_t

tic
for i=1:length(maxBins)
    for rep=1:replicates
        SPD_t(:,rep) = [metricPD(MANU_phenotypes{i,rep},'hamming'); metricPD(AUTO_phenotypes{i,rep},'hamming'); metricPD(AUTO_phenotypes5{i,rep},'hamming'); metricPD(AUTO_phenotypes10{i,rep},'hamming')];
        disp(['PDp: ' num2str(toc)]);
    end
    if replicates==1
        PDp(:,i) = (SPD_t');
    else
        PDp(:,i) = median(SPD_t');
    end
end
clear SPD_t

%% Fitness
clear fit;
for i=1:length(maxBins)
    fit(:,i) = [median([MANU_fitness{i,:}]); median([AUTO_fitness{i,:}]); median([AUTO_fitness5{i,:}]); median([AUTO_fitness10{i,:}])];
end

%% Analysis diversity over runs
manuMapGens = pStd.nGens:pStd.nGens:pStd.nGens*poemCfg.numIterations;
for i=1:length(maxBins)
    disp(['Max Bins: ' int2str(maxBins(i))]);
    for rep=1:replicates
        disp(['rep: ' int2str(rep)]);
        for t=1:length(stats{i,rep}.maps)
            disp(['t: ' int2str(t)]);
            
            genes = allMaps{i,rep}{manuMapGens(t)}.genes;
            [~,poly,~] = d.fitfun(genes);
            phenos = getPhenotypeBoolean(poly,d.resolution);
            PD_Iterations(1,i,t,rep) = metricPD(phenos,'hamming');
            fitIterMed(1,i,t,rep) = median(allMaps{i,rep}{manuMapGens(t)}.fitness);
            fitIter10(1,i,t,rep) = prctile(allMaps{i,rep}{manuMapGens(t)}.fitness,10);
            fitIter90(1,i,t,rep) = prctile(allMaps{i,rep}{manuMapGens(t)}.fitness,90);

            genes = stats{i,rep}.maps{t}.genes;
            [~,poly,~] = d.fitfun(genes);
            phenos = getPhenotypeBoolean(poly,d.resolution);
            PD_Iterations(2,i,t,rep) = metricPD(phenos,'hamming');
            fitIterMed(2,i,t,rep) = median(stats{i,rep}.maps{t}.fitness);
            fitIter10(2,i,t,rep) = prctile(stats{i,rep}.maps{t}.fitness,10);
            fitIter90(2,i,t,rep) = prctile(stats{i,rep}.maps{t}.fitness,90);

            genes = stats5{i,rep}.maps{t}.genes;
            [~,poly,~] = d.fitfun(genes);
            phenos = getPhenotypeBoolean(poly,d.resolution);
            PD_Iterations(3,i,t,rep) = metricPD(phenos,'hamming');
            fitIterMed(3,i,t,rep) = median(stats5{i,rep}.maps{t}.fitness);
            fitIter10(3,i,t,rep) = prctile(stats5{i,rep}.maps{t}.fitness,10);
            fitIter90(3,i,t,rep) = prctile(stats5{i,rep}.maps{t}.fitness,90);

            genes = stats10{i,rep}.maps{t}.genes;
            [~,poly,~] = d.fitfun(genes);
            phenos = getPhenotypeBoolean(poly,d.resolution);
            PD_Iterations(4,i,t,rep) = metricPD(phenos,'hamming');
            fitIterMed(4,i,t,rep) = median(stats10{i,rep}.maps{t}.fitness);
            fitIter10(4,i,t,rep) = prctile(stats10{i,rep}.maps{t}.fitness,10);
            fitIter90(4,i,t,rep) = prctile(stats10{i,rep}.maps{t}.fitness,90);

        end
    end
end


%% Visualization
%showPhenotype(MANU_genomes{r,rep},d)
%showPhenotype(AUTO_genomes{r,rep},d)

set(0,'DefaultFigureWindowStyle','default')
clear figs;
lw = 2;

% Mean Distance NN
figs(1) = figure(1);hold off;
semilogx(repmat(maxBins,4,1)',MODNN', 'LineWidth', lw);
grid on;% legend('Manual', 'Auto2', 'Auto5', 'Auto10' ,'Location','NorthWest'); title('Sum Genetic Distance to Nearest Neighbour');
xlabel('Max. # Bins');
axis([0 maxBins(end) 0 2000]);
ax = gca; ax.XTick = maxBins;

figs(2) = figure(2);hold off;
semilogx(repmat(maxBins,4,1)',MODNNp', 'LineWidth', lw);
grid on; %legend('Manual', 'Auto2', 'Auto5', 'Auto10' ,'Location','NorthWest'); title('Sum Phenotypic Distance to Nearest Neighbour');
xlabel('Max. # Bins');
axis([0 maxBins(end) 0 15]);
ax = gca; ax.XTick = maxBins;

% Solow-Polasky Diversity
figs(3) = figure(3);hold off;
semilogx(repmat(maxBins,4,1)',SPD', 'LineWidth', lw);
grid on; %legend('Manual', 'Auto2', 'Auto5', 'Auto10' ,'Location','NorthWest'); title('Genetic Solow-Polasky Diversity');
xlabel('Max. # Bins');
axis([0 maxBins(end) 0 350]);
ax = gca; ax.XTick = maxBins;

figs(4) = figure(4);hold off;
semilogx(repmat(maxBins,4,1)',SPDp', 'LineWidth', lw);
grid on; %legend('Manual', 'Auto2', 'Auto5', 'Auto10' ,'Location','NorthWest'); title('Phenotypic Solow-Polasky Diversity');
xlabel('Max. # Bins');
axis([0 maxBins(end) 0 300]);
ax = gca; ax.XTick = maxBins;

% Pure Diversity
figs(5) = figure(5);hold off;
semilogx(repmat(maxBins,4,1)',PD', 'LineWidth', lw);
grid on; %legend('Manual', 'Auto2', 'Auto5', 'Auto10' ,'Location','NorthWest'); title('Genetic Pure Diversity');
xlabel('Max. # Bins');
axis([0 maxBins(end) 0 15e13]);
ax = gca; ax.XTick = maxBins;

figs(6) = figure(6);hold off;
semilogx(repmat(maxBins,4,1)',PDp', 'LineWidth', lw);
grid on; %legend('Manual', 'Auto2', 'Auto5', 'Auto10' ,'Location','NorthWest'); 
title('Phenotypic Pure Diversity');
xlabel('Max. # Bins');
axis([0 maxBins(end) 0 15]);
ax = gca; ax.XTick = maxBins;


% Fitness
figs(end+1) = figure;hold off;
semilogx(repmat(maxBins,4,1)',fit', 'LineWidth', lw);
grid on; %legend('Manual', 'Auto2', 'Auto5', 'Auto10' ,'Location','NorthWest'); 
title('Avg. Fitness');
xlabel('Max. # Bins');
axis([maxBins(1) maxBins(end) 0.9 1.05]);
ax = gca; ax.XTick = maxBins;

save_figures(figs,'.','auto',11,[4 3])

%%

%% Fitness
set(0,'DefaultFigureWindowStyle','docked')

rep = 3
for r=1:length(maxBins)
    figs(end+1) = figure;hold off;
    plot(MANUmap{r,rep}.stats.fitnessMean);
    hold on;
    plot(AUTOmap{r,rep}.stats.fitnessMean);
    title(['Mean Fitness Gen']);
    xlabel('Gens'); ylabel('Mean Fitness (%)');
    legend('Manual','Auto 2','Auto 5','Auto 10','Location','SouthEast');
    ax = gca;
    ax.YLim = [0.9 1];
    grid on;
end

%%
rep=1;
for i=1:length(maxBins)
    
    showPhenotype(MANU_genomes{i,rep},d);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
    showPhenotype(AUTO_genomes{i,rep},d);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
    showPhenotype(AUTO_genomes5{i,rep},d);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
    showPhenotype(AUTO_genomes10{i,rep},d);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
end


%%
rep = 2;
for i=4:length(maxBins)
    nShapes = maxBins(i);
    xPos = 0:ceil(sqrt(nShapes))-1; [X,Y] = ndgrid(xPos,xPos);
    placement = [X(:) Y(:)]; placement = placement(1:nShapes,:);
    
    genomes = MANU_genomes{i,rep};
    [fitSort,idSort] = sort(MANU_fitness{i,rep},'ascend');
    genomes = genomes(idSort,:);
    showPhenotype(genomes,d,[],placement);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;

    genomes = AUTO_genomes{i,rep};
    [fitSort,idSort] = sort(AUTO_fitness{i,rep},'ascend');
    genomes = genomes(idSort,:);
    showPhenotype(genomes,d,[],placement);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;

    genomes = AUTO_genomes5{i,rep};
    [fitSort,idSort] = sort(AUTO_fitness5{i,rep},'ascend');
    genomes = genomes(idSort,:);
    showPhenotype(genomes,d,[],placement);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;

    genomes = AUTO_genomes10{i,rep};
    [fitSort,idSort] = sort(AUTO_fitness10{i,rep},'ascend');
    genomes = genomes(idSort,:);
    showPhenotype(genomes,d,[],placement);
    ax = gca; ax.XTick = []; ax.YTick = [];
    axis tight;
    
end


%% Diversity over runs
set(0,'DefaultFigureWindowStyle','docked')
manuMapGens = pStd.nGens:pStd.nGens:pStd.nGens*poemCfg.numIterations;


% PD_Iterations(1,i,t,rep)
PD_IterationsMean = mean(PD_Iterations,4);
figs(1) = figure;
plot( (squeeze(PD_IterationsMean(1,:,:)))' )
axis([1 2 0 20]);grid on;ylabel('PD Pheno');title('Manu');
xlabel('# Gens');ax = gca;ax.XTick = 1:poemCfg.numIterations; ax.XTickLabel = manuMapGens;
l = legend('50','100','200','400');l.Title.String = '# Bins';

figs(end+1) = figure;
plot( (squeeze(PD_IterationsMean(2,:,:)))' )
axis([1 2 0 20]);grid on;ylabel('PD Pheno');title('Auto 2');
xlabel('# Gens');ax = gca;ax.XTick = 1:poemCfg.numIterations; ax.XTickLabel = manuMapGens;
l = legend('50','100','200','400');l.Title.String = '# Bins';

figs(end+1) = figure;
plot( (squeeze(PD_IterationsMean(3,:,:)))' )
axis([1 2 0 20]);grid on;ylabel('PD Pheno');title('Auto 5');
xlabel('# Gens');ax = gca;ax.XTick = 1:poemCfg.numIterations; ax.XTickLabel = manuMapGens;
l = legend('50','100','200','400');l.Title.String = '# Bins';

figs(end+1) = figure;
plot( (squeeze(PD_IterationsMean(4,:,:)))' )
axis([1 2 0 20]);grid on;ylabel('PD Pheno');title('Auto 10');
xlabel('# Gens');ax = gca;ax.XTick = 1:poemCfg.numIterations; ax.XTickLabel = manuMapGens;
l = legend('50','100','200','400');l.Title.String = '# Bins';

% fitIterMed fitIter10 fitIter90
figs(end+1) = figure;hold off;
plot( median(squeeze(median(fitIterMed(1,:,:,:),4))),'k-'); hold on;
plot( median(squeeze(median(fitIter10(1,:,:,:),4))),'k--' );
plot( median(squeeze(median(fitIter90(1,:,:,:),4))),'k--' );
axis([1 2 0.8 1]);grid on;ylabel('Fitness');title('Manu');
legend('Median', '10%', '90%', 'Location', 'SouthEast');
xlabel('# Gens');ax = gca;ax.XTick = 1:poemCfg.numIterations; ax.XTickLabel = manuMapGens;

figs(end+1) = figure;hold off;
plot( median(squeeze(median(fitIterMed(2,:,:,:),4))),'k-' ); hold on;
plot( median(squeeze(median(fitIter10(2,:,:,:),4))),'k--' );
plot( median(squeeze(median(fitIter90(2,:,:,:),4))),'k--' );
axis([1 2 0.8 1]);grid on;ylabel('Fitness');title('Auto 2');
legend('Median', '10%', '90%', 'Location', 'SouthEast');
xlabel('# Gens');ax = gca;ax.XTick = 1:poemCfg.numIterations; ax.XTickLabel = manuMapGens;

figs(end+1) = figure;hold off;
plot( median(squeeze(median(fitIterMed(3,:,:,:),4))),'k-' ); hold on;
plot( median(squeeze(median(fitIter10(3,:,:,:),4))),'k--' );
plot( median(squeeze(median(fitIter90(3,:,:,:),4))),'k--' );
axis([1 2 0.8 1]);grid on;ylabel('Fitness');title('Auto 5');
legend('Median', '10%', '90%', 'Location', 'SouthEast');
xlabel('# Gens');ax = gca;ax.XTick = 1:poemCfg.numIterations; ax.XTickLabel = manuMapGens;

figs(end+1) = figure;hold off;
plot( median(squeeze(median(fitIterMed(4,:,:,:),4))),'k-' ); hold on;
plot( median(squeeze(median(fitIter10(4,:,:,:),4))),'k--' );
plot( median(squeeze(median(fitIter90(4,:,:,:),4))),'k--' );
axis([1 2 0.8 1]);grid on;ylabel('Fitness');title('Auto 10');
legend('Median', '10%', '90%', 'Location', 'SouthEast');
xlabel('# Gens');ax = gca;ax.XTick = 1:poemCfg.numIterations; ax.XTickLabel = manuMapGens;


%save_figures(figs,'.','divFitManuAuto',14,[5 5])

%%
model = configs{r,rep}{end}.model;
%sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
%samples = sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:);
samples = lhsdesign(100,2);

% Sample latent space equidistantly
minLat = min(latent);
maxLat = max(latent);

N = 20;
x = linspace(minLat(1), maxLat(1), N);
y = linspace(minLat(2), maxLat(2), N);

[X,Y] = ndgrid(x,y);

samplePts = [X(:),Y(:)];

returnImage = sampleAE(samplePts',poemCfg.model.decoderNet);
returnImage = squeeze(returnImage);
returnImage = permute(returnImage,[3,1,2]);

clear out;
for i=1:size(returnImage,1)
    out{i} = squeeze(returnImage(i,:,:));
end

figs(5) = figure(5);hold off;
out = imtile(out,'GridSize', [N N]);
imagesc(out);



