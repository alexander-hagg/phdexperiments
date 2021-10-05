clear;clc;
%%
DOF = 16;DOMAIN = 'catmullRom';ALGORITHM = 'grid';LATENTMODEL = 'VAE';
addpath(genpath('/home/alex/IPPM')); rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN])); rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]); rmpath(genpath('latentmodels')); addpath(genpath(['latentmodels/' LATENTMODEL]));
d = domain(DOF);
p = defaultParamSet;
d.resolution = 128;

m = cfgLatentModel('data/workdir',d.resolution);
poemCfg = poemParamSet(p,m);
poemCfg.numInitSamples                = p.numInitSamples;
poemCfg.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);
poemCfg.numIterations = 1;

%fname = ['catmullrom_asymmetry']; FITNESSFUNCTION = 'pointAsymmetry';  
fname = ['catmullrom_symmetry']; FITNESSFUNCTION = 'pointSymmetry';  

rmpath(genpath('domain/catmullRom/fitnessFunctions')); addpath(genpath(['domain/catmullRom/fitnessFunctions/' FITNESSFUNCTION]));

%% Experiment
% User Interaction Demonstration
p.numInitSamples = 64;

% Initialize initial solution set
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,polygons] = fitfun(initSamples,d);

[map{1}, config{1}, stats{1}] = poem(initSamples,polygons,fitness,poemCfg,d,2);
close all; clear figs;save([fname 'step1.mat']);
disp('Optimization Step 1 Done');


%%
load([fname 'step1.mat']);
[genes,fitness,features,bins] = extractMap(map{1});

%FITNESSFUNCTION = 'userANDasymmetry'; 
%selectionIDs = [507,531,534,557,558,559,560,561,585,586]; % asymmetry

FITNESSFUNCTION = 'userANDsymmetry'; 
%selectionIDs = [176,201,205,245]; % symmetry
selectionIDs = [598, 600, 627,  650, 741, 744, 746, 747, 748, 749]; % symmetry

rmpath(genpath('domain/catmullRom/fitnessFunctions')); addpath(genpath(['domain/catmullRom/fitnessFunctions/' FITNESSFUNCTION]));
d.userModel = stats{1}.models{end};



figs(2) = figure(4);
dSel = d;
dSel.phenoDistMult = 3;
showPhenotype(genes(selectionIDs,:),dSel,figs(2));title('Selected Shapes');
axis equal;



figs(1) = figure(5);
cmap = [0 0 0; 0 0 1];
colors = repmat(cmap(1,:),size(genes,1),1);
colors(selectionIDs,:) = repmat(cmap(2,:),numel(selectionIDs),1);
showPhenotype(genes,d,figs(1),bins,colors);


phenotypes = d.getPhenotype(genes);
features = predictFeatures(phenotypes,d.userModel);
d.selectedShapes = features(selectionIDs,:);
d.deselectedShapes = features; d.deselectedShapes(selectionIDs,:) = [];


%% Initialize new initial solution set

% Create mutated versions of selection
initMutationSigma = 0.1;

% Inject selected solutions
newSamples = genes(selectionIDs,:);
nNewPerSelected = ceil(p.numInitSamples./length(selectionIDs));
for i=1:length(selectionIDs)
    newSampleMutations = initMutationSigma * randn(nNewPerSelected,d.dof);
    newSamples = [newSamples; genes(selectionIDs(i),:) + newSampleMutations];
end

figs(3) = figure(6);
d.phenoDistMult = 10;
showPhenotype(newSamples,d,figs(3));title('Injected Shapes');


%initMutationSigma = 0.1;
% Inject selected solutions
%newSamples = genes(selectionIDs,:);
%nNewPerSelected = ceil(p.numInitSamples./length(selectionIDs));
%for i=1:length(selectionIDs)
%    newSampleMutations = initMutationSigma * randn(nNewPerSelected,d.dof)
%    newSamples = [newSamples; genes(selectionIDs(i),:) + newSampleMutations];
%end

%figs(99) = figure(99);
%d.phenoDistMult = 10;
%showPhenotype(newSamples,d,figs(99));title('Injected Shapes');

% 
% %% No constraint
% d.selectPenalty = 'none'; 
% [newSamplesfitness,newSamplespolygons] = fitfun(newSamples,d); % Recalculate fitness and phenotypes
% 
% [map{2}, config{2}, stats{2}] = poem(newSamples,newSamplespolygons,newSamplesfitness,poemCfg,d,2);
% 
% close all; clear figs; save([fname 'step2a.mat']);
% disp('Optimization Step 2a Done');
% 
% 
% %% inverseDistance
% d.selectPenalty = 'inverseDistance'; 
% [newSamplesfitness,newSamplespolygons] = fitfun(newSamples,d); % Recalculate fitness and phenotypes
% 
% [map{3}, config{3}, stats{3}] = poem(newSamples,newSamplespolygons,newSamplesfitness,poemCfg,d,2);
% 
% close all; clear figs; save([fname 'step2b.mat']);
% disp('Optimization Step 2b Done');
% 
% %% relativeDistance
% d.selectPenalty = 'relativeDistance';
% [newSamplesfitness,newSamplespolygons] = fitfun(newSamples,d); % Recalculate fitness and phenotypes
% 
% [map{4}, config{4}, stats{4}] = poem(newSamples,newSamplespolygons,newSamplesfitness,poemCfg,d,2);
% 
% close all; clear figs; save([fname 'step2c.mat']);
% disp('Optimization Step 2c Done');

%% relativeDistanceOnlyPenalizeConstraintViolation
d.selectPenalty = 'relativeDistanceOnlyPenalizeConstraintViolation'; 
[newSamplesfitness,newSamplespolygons] = fitfun(newSamples,d); % Recalculate fitness and phenotypes

[map{5}, config{5}, stats{5}] = poem(newSamples,newSamplespolygons,newSamplesfitness,poemCfg,d,2);

close all; clear figs; save([fname 'step2d.mat']);
disp('Optimization Step 2d Done');

%%
load([fname 'step2d.mat']);
d.phenoDistMult = 20;

for i=2:5
    [genes,fitness,~,bins] = extractMap(map{i});
    d.selectPenalty = 'relativeDistanceOnlyPenalizeConstraintViolation'; 
    [~,~,fitness] = fitfun(genes,d);
    fitness = fitness(:,2);
    %selected = fitness>=0.5;
    %genes = genes(selected,:);fitness = fitness(selected);bins = bins(selected,:);
    figs(2+i) = figure(5+i);
    distances = pdist2(genes,newSamples(1:numel(selectionIDs),:));
    [row,col] = find(distances==0);
    cmap = redgreencmap(256+1,'Interpolation','linear'); cmap = flipud(cmap);
    fitBins = floor(256*fitness)+1;
    colors = cmap(fitBins,:);
    colors(row,:) = repmat([0 0 1],numel(row),1);
    faceAlpha = fitness;
    faceAlpha(faceAlpha>0.5) = 1;
    showPhenotype(genes,d,figs(2+i),bins,colors,faceAlpha);
end




%% Analyze distance to selection
clear fitness;
for i=2:5
    [genes,~,~,bins] = extractMap(map{i});
    d.selectPenalty = 'relativeDistanceOnlyPenalizeConstraintViolation'; 
    [~,~,fitness{i-1}] = fitfun(genes,d);
    %selected = fitness(:,2)>=0.5;
    %genes = genes(selected,:);fitness = fitness(selected,1);bins = bins(selected,:);
    figs(8) = figure(11);
    subplot(1,4,i-1);
    boxplot(fitness{i-1});ax = gca;ax.YAxis.Limits = [0 1];grid on;
    %legend('Symmetry','Selection');
end

%%
save_figures(figs,'.',[fname],12,[4 4])

%% TEST CATMULLROM DOMAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%
% genome = [1 0.1 1 1 1 0.1 1 1,zeros(1,8)];
%
% for i=1:100
%     disp(i)
%     genome(1:8) = 0.1 + 0.9*rand(8,1);
%     genome(9:16) = -0.25 + 0.5*rand(8,1)
%     [shape,ctlPts] = getPhenotypeFFD(genome,d.base);
%     figure(1);
%     hold off;
%     plot(shape{1}.Vertices(:,1),shape{1}.Vertices(:,2));
%     hold on;
%     scatter(ctlPts{1}(:,1),ctlPts{1}(:,2))
%     axis([-0.5 0.5 -0.5 0.5]);
%     drawnow;
%     pause(0.3);
% end
%
% %%
% for i=1:100
%     disp(i)
%     genome(1:8) = 0.1 + 0.9*rand(8,1);
%     genome(9:16) = -0.25 + 0.5*rand(8,1)
%     shape = getPhenotypeFFD(genome,d.base);
%     [~,booleanMap] = getPhenotypeBoolean(shape)
%     figure(1);
%     imagesc(booleanMap{1})
%     drawnow;
%     pause(0.1);
% end



