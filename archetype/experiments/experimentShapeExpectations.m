clear;clc;
DOF = 16;
DOMAIN = 'footprints';
ALGORITHM = 'voronoi';

addpath(genpath('/home/alex/archetype'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);

d = domain(DOF);
p = defaultParamSet;
d.fitfun = d.fitfunPointSymm; % Multimodal function

scale = 0.2:0.1:1

%% Run simulation
clear basicShapes;

experimentName = 'basicShapes';
for i=1:numel(scale)
    sc = scale(i);
    % Small to large circles
    basicShapes(i,:) = [sc.*(ones(1,d.dof/2)), zeros(1,d.dof/2)];
end
for i=1:numel(scale)
    sc = scale(i);
    % Star to circle
    basicShapes(numel(scale)+i,:) = [reshape([ones(1,d.dof/4);sc.*ones(1,d.dof/4)],1,d.dof/2), zeros(1,d.dof/2)];
end

d.nGPUs = 1;
[fitness,features,polyshapes,booleanMap,rawfeatures,allData] = fitnessLettuce(basicShapes,d)

save([experimentName '.mat']);
%% Holdout Crossvalidation of GP models
clear gpModel p prediction output input
p = paramsGP(d.dof);

%trainData = basicShapes(1:9,:);
%allFeatures = rawfeatures(1:9,:);
trainData = basicShapes(10:18,:);
allFeatures = rawfeatures(10:18,:);

for i=1:numel(scale)
    disp(['Run ' int2str(i) '/' int2str(numel(scale))]);
    input = trainData;
    input(i,:) = [];
    
    output = allFeatures(:,3);
    output(i) = [];
    [gpModel] = trainGP(input,output,p);
    pred = predictGP(gpModel, trainData(i,:));
    prediction(i,1,:) = pred;
    
    output = allFeatures(:,4);
    output(i) = [];
    [gpModel] = trainGP(input,output,p);
    pred = predictGP(gpModel, trainData(i,:));
    prediction(i,2,:) = pred;
    
end


mape(1) = mean(abs(prediction(:,1,1)-allFeatures(:,3))./allFeatures(:,3))*100;
mape(2) = mean(abs(prediction(:,2,1)-allFeatures(:,4))./allFeatures(:,4))*100;
disp(['Mean Absolute Percentage Errors'])
disp(['E: ' num2str(mape(1))])
disp(['Umax: ' num2str(mape(2))])

%% Visualization
clear figs;close all

fig(1) = figure(1);
yyaxis left
hold off;
plot(allFeatures(:,3),'-');
hold on;
plot(prediction(:,1,1),'--')
plot(prediction(:,1,1)+20*prediction(:,1,2),':')
plot(prediction(:,1,1)-20*prediction(:,1,2),':')
ylabel('Enstrophy');
%
yyaxis right
hold off;
plot(allFeatures(:,4),'-');
hold on;
plot(prediction(:,2,1),'--')
plot(prediction(:,2,1)+20*prediction(:,2,2),':')
plot(prediction(:,2,1)-20*prediction(:,2,2),':')
ylabel('uMax');
grid on;
ax = gca;
ax.XTick = [1:9];
title(experimentName);
ax.XAxis.Limits = [0.5 9.5];
legend('True','Prediction','+/- 20*Variance');

save_figures(fig, '.', experimentName, 12, [6 4]);

%% Parameterize with Holdout Crossvalidation of GP models
clear gpModel p prediction output input loglikelihoodEvolution pred mape
p = paramsGP(size(trainData,2));

% name = 'circles'; trainData = basicShapes(1:9,:); allFeatures = rawfeatures(1:9,:);
name = 'stars'; trainData = basicShapes(10:18,:); allFeatures = rawfeatures(10:18,:);

lengthscales = [-3 -2 -1 0 1 2 3];
noiseVar = [-2 -1 0 1 2];
l = 1;
for k=1:length(lengthscales)
    for j=1:length(noiseVar)
        
        for i=1:numel(scale)
            disp(['Run ' int2str(i) '/' int2str(numel(scale))]);
            input = trainData;
            input(i,:) = [];
            
            output = allFeatures(:,3);
            output(i) = [];
            p.hyp.mean = mean(output);
            
            p.hyp.cov(1) = lengthscales(k);
            p.hyp.cov(2) = noiseVar(j);
            
            disp(['Initial length scale log(' num2str(10^p.hyp.cov(1)) ')']);
            disp(['Initial noise var log(' num2str(10^p.hyp.cov(2)) ')']);
            [gpModel] = trainGP(input,output,p,false);
            
            pred = predictGP(gpModel, trainData(i,:));
            prediction(i,j,k,:) = pred;
            loglikelihood(i,j,k) = gpModel.loglikelihoodEvolution(end);
        end
        mape(j,k) = mean(abs(prediction(:,j,k,1)-allFeatures(:,3))./allFeatures(:,3))*100
    end
end


save([name '.mat']);

%%
clear fig; close all;
set(0,'DefaultFigureWindowStyle','default')

fig(1) = figure(1);
imagesc(squeeze(mean(loglikelihood,1)))
cb = colorbar; cb.Label.String = 'Log Likelihood';
caxis([-75 -60]);
ax = gca; 
ax.XTick = 1:length(lengthscales); ax.XTickLabel = lengthscales;
ax.YTick = 1:length(noiseVar); ax.YTickLabel = noiseVar;
xlabel('Length scale'); ylabel('Noise \sigma');

fig(2) = figure(2);
imagesc(log(mape))
cb = colorbar; cb.Label.String = 'log(MAPE)';
caxis([log(0) log(40)]);
ax = gca; 
ax.XTick = 1:length(lengthscales); ax.XTickLabel = lengthscales;
ax.YTick = 1:length(noiseVar); ax.YTickLabel = noiseVar;
xlabel('Length scale'); ylabel('Noise \sigma');





%%
save_figures(fig, '.', [name '_parameterization'], 12, [5 3]);


