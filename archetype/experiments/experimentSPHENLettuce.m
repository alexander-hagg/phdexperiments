%experimentSPHENLettuce - Run Lettuce CFD domain
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Apr 2020; Last revision: 11-Apr-2019
%
%------------- BEGIN CODE --------------

clear;clc;
DOF = 16;
workpath = ['/home/' getenv('USER') '/archetype/'];
addpath(genpath(workpath));
DOMAIN = 'footprints';  rmpath(genpath([workpath 'domain'])); addpath(genpath([workpath 'domain/' DOMAIN]));
QD = 'grid';            rmpath(genpath([workpath 'QD/grid'])); rmpath(genpath([workpath 'QD/voronoi'])); addpath(genpath([workpath 'QD/' QD]));
SAQD = 'sphen';          rmpath(genpath([workpath 'QD/sail'])); rmpath(genpath(['QD/sphen'])); addpath(genpath([workpath 'QD/' SAQD]));

nGPUs = str2num(getenv('numGPUs')); if isempty(nGPUs); nGPUs = 4; end
d = domain(DOF,'lettuce',nGPUs);
p = defaultParamSet;

p.infill = infillParamSet;
p.infill.nTotalSamples = 1024;
p.infill.nAdditionalSamples = 16;
p.infill.intermediateMaps = false;
p.infill.injectTrueValues = true;

sphenP = p; sphenP.nGens = 2^10;  sphenP.nChildren = 2^5; sphenP.resolution = 32;

d.fitfun = d.fitfunLettuce; % Multimodal function: point symmetry (center); for testing purposes
experimentName = 'LETTUCE'
numReplicates = 1;

%%
initSampleSet = 'initsamples64.mat';
if exist(initSampleSet,'file')
    load('initsamples64.mat');
    for i=1:4
        features(:,i) = (initSet.unnormalizedFeatures(:,i)-d.featureMin(i))./(d.featureMax(i)-d.featureMin(i));
    end
    features(features>1) = 1; features(features<0) = 0;
    initSet.features = features;
    initSet.fitness = (1./(1+initSet.features(:,4)))*2-1;
    initSet.features = initSet.features(:,d.selectedFeatures);
else    
    disp('Creating initial sample set');
    sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen'); sobPoint = 1;
    initSet.samples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
    [initSet.fitness,~,~,~,initSet.unnormalizedFeatures,~,initSet.features] = d.fitfun(initSet.samples);
    
    for i=1:4
        features(:,i) = (initSet.unnormalizedFeatures(:,i)-d.featureMin(i))./(d.featureMax(i)-d.featureMin(i));
    end
    features(features>1) = 1; features(features<0) = 0;
    initSet.features = features;
    initSet.fitness = (1./(1+initSet.features(:,4)))*2-1;
    save('initsamples64.mat','initSet');
    initSet.features = initSet.features(:,d.selectedFeatures);
end
%% ----------------------------------------------------------------------------------

for rep=1 : numReplicates
    %% SPHEN (Surrogate-Assisted Phenotypic Niching
    [mapSPHEN{rep},surrogateFitnessSPHEN{rep},surrogateFeaturesSPHEN{rep},allMapsSPHEN{rep},trueFilledSPHEN] = sphen(initSet,sphenP,d,1);    
    save([experimentName '_replicate_' int2str(rep) '.mat']);
end

%% Analysis
% set(0,'DefaultFigureWindowStyle','default')

% rep = 1;
% [medresmap{rep},errors(1,rep,:),filled(1,rep),medianFitness(1,rep)] = analyzeMaps(mapSPHEN{rep},d,p,false,16);
% [lowresmap{rep},errors(1,rep,:),filled(1,rep),medianFitness(1,rep)] = analyzeMaps(mapSPHEN{rep},d,p,false,3);


% fig(1) = figure(1);
% viewMap(medresmap{rep},d);
% caxis([0.88 0.92]);

% fig(2) = figure(2);
% viewMap(lowresmap{rep},d)
% caxis([0.88 0.92]);

% truefit = [0.9240    0.9084    0.8905    0.9116    0.9031    0.8923    0.9032    0.9005    0.8920];
% truefitmap = lowresmap{rep};
% truefitmap.fitness = reshape(truefit,3,3);
% fig(3) = figure(3);
% viewMap(truefitmap,d)
% caxis([0.88 0.92]);

% save_figures(fig, '.', ['LettuceSPHEN'], 12, [5 4]);


%%
% tRes = size(medresmap{rep}.fitness,1);
% x = 1:tRes;
% y = 1:tRes;
% [X,Y] = ndgrid(x,y);

% genes = reshape(medresmap{rep}.genes,[],d.dof);
%fitness = reshape(medresmap{rep}.fitness,[],1);
% features = reshape(medresmap{rep}.features,[],2);

% fig(4) = showPhenotype(genes,d,[],[X(:),Y(:)])
% axis equal;

% tRes = size(lowresmap{rep}.fitness,1);
% x = 1:tRes;
% y = 1:tRes;
% [X,Y] = ndgrid(x,y);

% genes = reshape(lowresmap{rep}.genes,[],d.dof);
% fitness = reshape(lowresmap{rep}.fitness,[],1);
% features = reshape(lowresmap{rep}.features,[],2);

% fig(5) = showPhenotype(genes,d,[],[X(:),Y(:)])
% axis equal;



%%
% fig(6) = showPhenotype(surrogateFitnessSPHEN{rep}.trainInput,d)
% axis equal;
%%
% in = surrogateFitnessSPHEN{rep}.trainInput;
% in = mapminmax(in',0,1)';
% %out = surrogateFitnessSPHEN{rep}.trainOutput;
% out = surrogateFeaturesSPHEN{rep}{1}.trainOutput;
% out = mapminmax(out',0,1)';
% len = [600];
% randSet = rand(1000,16);
% %figure(1);hold off;
% for i=1:length(len)
% nModel = trainGP(in(1:len(i),:),out(1:len(i),1),paramsGP(16));
% [prediction] = predictGP(nModel,randSet);
% mean(prediction(:,2))
 %std(prediction(:,2))
 %%subplot(length(len),1,i);
 %%boxplot(prediction(:,2));
 %end
%%
% figure(2);hold off;
%scatter(1:length(out),1./out-1);
%hold on;
%plot(movmean(1./out-1,100))
% scatter(1:length(out),out);
% hold on;
% plot(movmean(out,100))

%%
% save_figures(fig, '.', ['LettuceSPHEN'], 12, [5 4]);

%%
% workdir = '/home/alex/archetype/data/lettuce/';

% pgons = d.getPhenotype(genes);
% [~,images] = getPhenotypeBoolean(pgons);
 %
% for i=1:length(images)
%     mkdir([workdir int2str(i)]);
%     cd([workdir int2str(i)]);
%     imwrite(images{i},['building.png'])
 %    system('cp /home/alex/archetype/domain/footprints/exampleCallLettuce/building.py .');
%     system('cp /home/alex/archetype/domain/footprints/exampleCallLettuce/DMD_multiple.py .');
%     cd('..')
% end

%% True map
% mapSPHEN{rep}
