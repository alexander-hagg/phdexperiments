% EXPERIMENTDIVERSITY
%
% Chaos script
% 
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2020; Last revision: 03-Feb-2020
%
%------------- BEGIN CODE --------------

clear;clc;
init;

%% I) Multiobjective optimization with NSGA-II
numObjectives = evaluate_objective;
nsgaFun = @(x) evaluate_objective(x,d);
pc = 0.9; pm = 1/d.dof;
tic

[R, Rfit, Rrank, stats] = NSGAII(p.nChildren,p.nGens,pc,pm,p.mutSigma,nsgaFun,d.dof,d.ranges(:,1),d.ranges(:,2),initSamples);
% initSamples will be used for half the population (parents + children), so
% inside NSGAII I have changed the first rank selection to only select half
% as many parents as was provided in the initial sample set.
toc

select = Rrank==1;
%NSGA_genomes = R(select,:); %NSGA_fitness = Rfit(select,1)'; NSGA_features = [Rfit(select,2),1-Rfit(select,3)];
NSGA_genomes = R;
[NSGA_fitness,NSGA_polygons,~] = d.fitfun(NSGA_genomes);
NSGA_features = categorize(NSGA_polygons, d)';

NSGAmap = createMap(d, p);
[replaced, replacement, features]       = nicheCompete(NSGA_genomes, NSGA_fitness, NSGAmap, d, p, NSGA_features');
NSGAmap                                 = updateMap(replaced,replacement,NSGAmap,NSGA_fitness,NSGA_genomes,NSGA_features');

%% II) Illumination with QD (default)
%global fcnEvalCtrQD;
%fcnEvalCtrQD = 0;
p.selectProcedure = 'random'; % 'curiosity'
p.display.illu = false;

features = categorize(polygons, d);
map = createMap(d, p,p.extraMapValues);
[replaced, replacement] = nicheCompete(initSamples,fitness,map,d,p,features);
map = updateMap(replaced,replacement,map,fitness,initSamples,features);

tic
[QDmap, percImproved, percValid, allMaps, percFilled, fitnessMean, fitnessTotal] = illuminate(map,d.fitfun,p,d);
toc

%% III) Niching in Genetic Space

% Run python script and load results from disk (Domain is defined in Python script!)
numOuterIters = (p.featureResolution(1)*p.featureResolution(2));
numInnerEvals = ceil(p.nChildren*p.nGens/numOuterIters);
localSearchMethod = 'BFGS';
cd('MMO/basinhopping');
system(['python3 RLS.py -n ' int2str(numInnerEvals) ' -i ' int2str(numOuterIters) ' -l ' localSearchMethod ' -a ' num2str(axialMultiplier) ' -m ' num2str(radialMultiplier)]);
cd('../..');

%RLS_genomes = [csvread('MMO/basinhopping/optima.csv'); csvread('MMO/basinhopping/visited_points.csv')];
RLS_genomes = [csvread('MMO/basinhopping/optima.csv')];
[RLS_fitness,RLS_polygons] = d.fitfun(RLS_genomes);
RLS_features = categorize(RLS_polygons, d);
RLSmap = createMap(d, p);
[replaced, replacement, features]       = nicheCompete(RLS_genomes, RLS_fitness, RLSmap, d, p, RLS_features);
RLSmap                                 = updateMap(replaced,replacement,RLSmap,RLS_fitness',RLS_genomes,RLS_features);

%% Visualize
clear figs
set(0,'DefaultFigureWindowStyle','default')

viewMap(NSGAmap,d);caxis([0 1]);
title('MOO: NSGA-II in Morphological Feature Space');xlabel(d.featureLabels{1});ylabel(d.featureLabels{2});
figs(1) = gcf;
viewMap(QDmap,d);caxis([0 1]);
title('QD: MAP-Elites in Morphological Feature Space');xlabel(d.featureLabels{1});ylabel(d.featureLabels{2});
figs(3) = gcf;
viewMap(RLSmap,d);caxis([0 1]);
title('MMO: RLS in Morphological Feature Space');xlabel(d.featureLabels{1});ylabel(d.featureLabels{2});
figs(2) = gcf;

drawnow;


%% Visualize shapes

QDgenes = reshape(QDmap.genes,[],d.dof); QDgenes(all(isnan(QDgenes)'),:) = [];
allGenes = [NSGA_genomes;QDgenes;RLS_genomes];
ids = [ones(size(NSGA_genomes,1),1); 2*ones(size(QDgenes,1),1)];ids = [ids; 3*ones(size(RLS_genomes,1),1)];

polyshapes = d.getPhenotype(allGenes);
[flatPhenotypes,booleanMap] = getPhenotypeBoolean(polyshapes,64);

[coeffPhenotypic,scorePhenotypic,latentPhenotypic,tsquaredPhenotypic,explainedPhenotypic,muPhenotypic] = pca(flatPhenotypes);
[coeffGenetic,scoreGenetic,latentGenetic,tsquaredGenetic,explainedGenetic,muGenetic] = pca(allGenes);

[allFitness] = d.fitfun(allGenes);
%manifold = getAEConfig('data/workdir');manifold = manifold.train(booleanMap);simX = manifold.latent';
%%
cmap = [0 0 0; 1 0 0; 0 0 1];
spaces = {'','','Genetic','Phenotypic'};
for j=3:4
    simX = eval(['score' spaces{j} '(:,1:2)']);
    figs(j) = figure(j);hold off;
    
    if range(simX(:,2)) > range(simX(:,1))
        simX = fliplr(simX);
    end
    
    clear h
    skip = 0.1;
    minGrid = min(simX);maxGrid = max(simX);
    [xq,yq] = meshgrid(minGrid(1):skip:maxGrid(1), minGrid(2):skip:maxGrid(2));
    vq = griddata(simX(:,1),simX(:,2),allFitness-1,xq,yq,'natural');
    sf = surf(xq,yq,vq);
    sf.EdgeAlpha = 0;
    view(0,90);
    
    hold on;
    for i=1:3
        h(i) = scatter(simX(ids==i,1),simX(ids==i,2),16,cmap(ids(ids==i),:),'filled');
    end
    cb = colorbar; cb.Label.String = 'Fitness';
    title([spaces{j} ' Similarity (PCA, ' num2str(eval(['sum(explained' spaces{j} '(1:2))'])) ' % explained)']);
    l = legend(h,'NSGA-II (MOO)','MAP-Elites (QD)','RLS (MMO)','Location','NorthWest');
    ax = gca; ax.XTick = [];ax.YTick = [];
    axis equal; axis tight
end
%% Show polygons in phenotypic similarity space coordinates
simX = eval(['score' spaces{4} '(:,1:2)']);
    
placement = 2*simX;
cmap = parula(numel(unique(ids)));

lbls = {'NSGA-II','MAP-Elites','RLS'};
figs(102) = figure(102);hold off;
for i=1:max(unique(ids))
    %subplot(3,1,i)
    ax = gca;
    
    %
    showPhenotype(allGenes(ids==i,:),d,figs(102),placement(ids==i,:),cmap(ids(ids==i),:));
    hold on;
    axis equal
    axis([min(placement(:,1)) max(placement(:,1)) min(placement(:,2)) max(placement(:,2))]);
    title(lbls{i});
    
end
%%
lbls = {'NSGA-II','MAP-Elites','RLS'};

figs(101) = figure(101);hold off;
for i=1:max(unique(ids))
    subplot(3,1,i)
    histogram(allFitness(ids==i),[0:0.05:1]);
    title(lbls{i});
    xlabel('Fitness');
    ylabel('Count');
    axis([-1.1 1.1 0 p.numInitSamples/2]);
    grid on;
end

%title('NSGA-II, RLS and MAP-Elites phenotypes in similarity space')
%legend([ax.Children(1) ax.Children(end/2) ax.Children(end)],'NSGA-II (MOO)','MAP-Elites (QD)','RLS (MMO)')

%% Show polygons in archive coordinates

[fitnessAll,polygonsAll] = d.fitfun(allGenes);
featuresAll = categorize(polygonsAll, d);
mapAll = createMap(d, p,p.extraMapValues);
[replaced, replacement] = nicheCompete(allGenes,fitnessAll,mapAll,d,p,featuresAll);
mapAll = updateMap(replaced,replacement,mapAll,fitnessAll',allGenes,featuresAll);


cmap = parula(numel(unique(ids)));

figHandle = figure(99);hold off;
placement = 2*simX;
%
showPhenotype(allGenes,d,figHandle,placement,cmap(ids,:));
axis([min(placement(:,1)) max(placement(:,1)) min(placement(:,2)) max(placement(:,2))]);
axis equal
axis tight


%%

save_figures([figs(1) figs(2) figs(3)],'.','3algs',14,[8 8])
save_figures([figs(4) figs(5)],'.','3algs',14,[8 4])
save_figures([figs(100) figs(101)],'.','3algs',14,[6 9])


%% Parameter ranges of optimal set
%NSGA_genomes NSGA_fitness RLS_genomes RLS_fitness genes = reshape(QDmap.genes,[],d.dof)
QD_genomes = reshape(QDmap.genes,[],d.dof);

ranges = [range(NSGA_genomes);range(QD_genomes);range(RLS_genomes)];
mins = [min(NSGA_genomes);min(QD_genomes);min(RLS_genomes)];
maxs = [max(NSGA_genomes);max(QD_genomes);max(RLS_genomes)];


figure(1);hold off;
for i=1:3
    pgon = polyshape([1:16 16:-1:1 1],[mins(i,:) fliplr(maxs(i,:)) mins(i,1)])
    plot(pgon,'FaceColor',cmap(i,:),'FaceAlpha',0.1);
    hold on;
end

legend('NSGA-II','QD','RLS');


%% Convex hull
% 'Ts' print out memory statistics
% 'T1' races the overall execution of the program

% Genomes
tic
[~,vol] = convhulln(NSGA_genomes(1:100,:),{'Qt','Qx','Ts'})
toc

% Phenotypes

% Monte Carlo f?r Volumbestimmung

%% Sweeps



%% END OF CODE