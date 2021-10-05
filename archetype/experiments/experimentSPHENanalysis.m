clear;clc;
DOF = 16;
DOMAIN = 'footprints';
QD = 'grid';
SAQD = 'sphen';

addpath(genpath('/home/alex/archetype'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath(genpath('QD/grid')); rmpath(genpath('QD/voronoi')); addpath(genpath(['QD/' QD]));
rmpath(genpath('QD/sail')); rmpath(genpath('QD/sphen')); addpath(genpath(['QD/' SAQD]));

d = domain(DOF);
p = defaultParamSet;
%p.intermediateMaps = true;

p.infill = infillParamSet;
p.infill.nTotalSamples = 1024;
p.infill.nAdditionalSamples = 16;%p.nChildren;

p.numInitSamples = 16;

sphenP = p;
sphenP.nGens = 2^9;
sphenP.nChildren = 2^6;

d.fitfun = d.fitfunPointSymm; % Multimodal function: point symmetry (center); for testing purposes
experimentName = 'SPHEN_analysis_pointSymmetry'
numReplicates = 3;
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');

%% ----------------------------------------------------------------------------------
set(0,'DefaultFigureWindowStyle','docked')

for rep=1 : numReplicates
    % Create initial sample set
    sobPoint = (rep-1)*p.numInitSamples + 1;
    initSet{rep}.genes = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
    [initSet{rep}.fitness,initSet{rep}.features] = d.fitfun(initSet{rep}.genes);

    %% SPHEN (Surrogate-Assisted Phenotypic Niching
    sphenP.infill.acqFcn               = 'PureUCB'; %FeatureUnCertainty FeatureCertainty PureUCB
    [mapSPHEN{rep},surrogateFitnessSPHEN{rep},surrogateFeaturesSPHEN{rep},allMapsSPHEN{rep},trueFilledSPHEN(rep)] = sail(initSet{rep},sphenP,d,1);
    %%
    sphenP.infill.acqFcn               = 'FeatureCertainty'; %FeatureUnCertainty FeatureCertainty PureUCB
    [mapSPHEN2{rep},surrogateFitnessSPHEN2{rep},surrogateFeaturesSPHEN2{rep},allMapsSPHEN2{rep},trueFilledSPHEN2(rep)] = sail(initSet{rep},sphenP,d,10);
    %%
    sphenP.infill.acqFcn               = 'FeatureUnCertainty'; %FeatureUnCertainty FeatureCertainty PureUCB
    [mapSPHEN3{rep},surrogateFitnessSPHEN3{rep},surrogateFeaturesSPHEN3{rep},allMapsSPHEN3{rep},trueFilledSPHEN3(rep)] = sail(initSet{rep},sphenP,d,20);
        
end

save([experimentName '_7.mat']);

%% Analysis
% SPHEN without feature acq
p.resolution = 8;
for rep=1:numReplicates
    [truemap{1,rep},errors(1,rep,:),filled(1,rep),medianFitness(1,rep)] = analyzeMaps(mapSPHEN{rep},d,p);
    [truemap{2,rep},errors(2,rep,:),filled(2,rep),medianFitness(2,rep)] = analyzeMaps(mapSPHEN2{rep},d,p);
    [truemap{3,rep},errors(3,rep,:),filled(3,rep),medianFitness(3,rep)] = analyzeMaps(mapSPHEN3{rep},d,p);
    
    for a=1:length(allMapsSPHEN{rep})
        disp(a)
        [alltruemap{1,a,rep},allerrors(1,a,rep,:),allfilled(1,a,rep),allmedianFitness(1,a,rep)] = analyzeMaps(allMapsSPHEN{rep}{a},d,p);
        [alltruemap{2,a,rep},allerrors(2,a,rep,:),allfilled(2,a,rep),allmedianFitness(2,a,rep)] = analyzeMaps(allMapsSPHEN2{rep}{a},d,p);
        [alltruemap{3,a,rep},allerrors(3,a,rep,:),allfilled(3,a,rep),allmedianFitness(3,a,rep)] = analyzeMaps(allMapsSPHEN3{rep}{a},d,p);
    end
end

%%
fig(1) = figure(1);hold off;
meanallfilled = mean(allfilled,3);
for i=1:3; plot(meanallfilled(i,:)); hold on; end
legend('PureUCB','FeatureCertainty','FeatureUnCertainty','Location','SouthEast');
title('Filled %');grid on;

fig(2) = figure(2);
titles = {'Fitness (Symmetry)','Area','Perimeter'};
meanallerrors = squeeze(mean(allerrors,3));
for i=1:3
    subplot(3,1,i);
    for j=1:3
        plot(meanallerrors(j,:,i)); 
        hold on; 
    end
    ax = gca; ax.YAxis.Limits = [0 0.12];
    title(titles{i});grid on;
    if i==1; legend('PureUCB','FeatureCertainty','FeatureUnCertainty'); end    
end


fig(3) = figure(3);hold off;
meanallmedianFitness = mean(allmedianFitness,3);
for i=1:3; plot(meanallmedianFitness(i,:)); hold on; end
legend('PureUCB','FeatureCertainty','FeatureUnCertainty','Location','SouthEast');
title('Median Fitness');grid on;


save_figures(fig, '.', [experimentName 'fitANDfill'], 12, [5 4]);