
clear;clc;
DOF = 16;DOMAIN = 'npoly_ffd_imagePhenotype';
addpath(genpath('/home/alex/poem'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));

ALGORITHM = 'grid';
rmpath(genpath('QD')); addpath(genpath(['QD/' ALGORITHM])); addpath(genpath(['QD/poem']));

d = domain(DOF);
p = poemParamSet;
%d.fitfun = d.fitfunAreaCirc; % Unimodal function
d.fitfun = d.fitfunPointSymm; % Multimodal function

%m = getAEConfig('data/workdir');
m = getTSNEConfig(DOF);

sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,phenotypes] = d.fitfun(initSamples);

selectionMethods = {'random','maxuncertainty'};
%%
numReplicates = 1;

for sel=1:length(selectionMethods)
    p.selectionMethod = selectionMethods{sel};
    for repl=1:numReplicates
        disp(['Selection method: ' p.selectionMethod]);
        disp(['Repetition ' int2str(repl) '/' int2str(numReplicates)]);
        profile on
        [map{sel,repl},phenoModel{sel,repl},stats{sel,repl}] = poem(initSamples,fitness,phenotypes,p,d,m);
        profile off
        profile viewer
    end
end

%% TODO Show movement of final map in tSNE
%% maybe.. TODO novelty selection (furthest from training examples
%% maybe.. TODO novelty plus curiosity

%% Show polygon maps 

cmap = parula(length(selectionMethods));
for i=1:length(selectionMethods)
    figHandle = figure(i);hold off;
    for j=1
        features = reshape(stats{i,j}.maps{20}.features,[],2);
        nanIDs = all(isnan(features)');
        features(nanIDs,:) = [];
        placement = 64*features;
        genes = reshape(stats{i,j}.maps{20}.genes,[],d.dof);
        genes(nanIDs,:) = [];
        showPhenotype(genes,d,figHandle,placement);  
        drawnow;
        hold on;
    end
end

%%
for i=1:length(selectionMethods)
    figHandle = figure(i);
    features = reshape(map{i}.features,[],2);
    placement = 64*features;
    genes = reshape(map{i}.genes,[],d.dof);
    genes(all(isnan(genes)'),:) = [];
    showPhenotype(genes,d,figHandle,placement);    
    title(selectionMethods{i});
    drawnow;
end
%%
cmap = parula(length(selectionMethods));
figure(5);hold off;
for i=1:length(selectionMethods)
    for j=1:5
        h(i,j) = plot(stats{i,j}.error.median,'Color',cmap(i,:))
        hold on;
    end
end
legend(h(:,1),selectionMethods);
title('Median Error');
%%

figure(6);hold off;
for i=1:length(selectionMethods)
    for j=1:5
        h(i,j) = plot(stats{i,j}.fitness.median,'Color',cmap(i,:))
        hold on;
    end
end
legend(h(:,1),selectionMethods);
axis([1 length(stats{i,j}.fitness.mean) 0 1]);
title('Median Fitness');

%%
figure(7);hold off;
for i=1:length(selectionMethods)
    for j=1:5
        h(i,j) = plot(stats{i,j}.elites.number,'Color',cmap(i,:))
        hold on;
    end
end
legend(h(:,1),selectionMethods);
title('Number of Solutions');

%%
clear h
figure(8);hold off;
for i=1:length(selectionMethods)
    %fitTotal = zeros(1,length(stats{i,j}.fitness.mean));
    fitTotal = [];
    for j=1:5
        fitTotal = [fitTotal;stats{i,j}.fitness.total];
    end
    h(i) = plot(median(fitTotal),'Color',cmap(i,:));     
    hold on;
end
legend(h,selectionMethods);
axis([1 length(stats{i,j}.fitness.mean) 0 300]);
title('Total Fitness');
%%

for i=1:length(selectionMethods)
    figHandle = figure(10+i);
    ax=gca;
    viewMap(map{i,1},d,ax);
    title(selectionMethods{i});
    drawnow;
end


%% TSNE
allGenes = [];
ids = [];
for i=1:length(selectionMethods)
    for j=1:5
        genes = reshape(map{i,j}.genes,[],d.dof);
        genes(all(isnan(genes)'),:) = [];
        allGenes = [allGenes;genes];
        ids = [ids, repmat([i],1,size(genes,1))];
    end
end

[~,~,~,phenotypes] = d.getPhenotype(allGenes);

simX = getSimSpace(phenotypes);
%%
scatter(simX(:,1),simX(:,2),32,cmap(ids,:),'filled')
%%
cmap = parula(length(selectionMethods));
figHandle = figure(99);hold off;
placement = 4*simX;
showPhenotype(allGenes,d,figHandle,placement,cmap(ids,:));  
    

%% END OF CODE