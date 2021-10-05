% testLettuce - test Lettuce solver
%
% Please include the following references in any publication using this code. For Bibtex please see the end of this file.
%
% Hagg, A., Asteroth, A. and B??ck, T., 2018, September. Prototype discovery using quality-diversity. In International Conference on Parallel Problem Solving from Nature (pp. 500-511). Springer, Cham.
% Hagg, A., Asteroth, A. and B??ck, T., 2019, July. Modeling user selection in quality diversity. In Proceedings of the Genetic and Evolutionary Computation Conference (pp. 116-124). ACM.
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2020; Last revision: 07-Jan-2020

%------------- BEGIN CODE --------------

clear;clc;
DOF = 16;
DOMAIN = 'footprints';
ALGORITHM = 'grid';

addpath(genpath('.'));
rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));
d = domain(DOF);

rmpath('QD/grid'); rmpath('QD/voronoi'); addpath(['QD/' ALGORITHM]);
p = defaultParamSet(4);

system('conda init bash')
system('conda activate lettuce');
    
%% Test
numInitSamples = 100;
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+numInitSamples)-1,:)+d.ranges(:,1)';


%%
clear initSamples;
scale = 0.2:0.05:1
%scale = [1.0 0.7 0.5 0.3]
d.resolution = 128;

for i=1:numel(scale)
    sc = scale(i);
    % Small to large circles
    %initSamples(i,:) = [sc.*(ones(1,d.dof/2)), zeros(1,d.dof/2)];
    % Star to circle
    initSamples(i,:) = [reshape([ones(1,d.dof/4);sc.*ones(1,d.dof/4)],1,d.dof/2), zeros(1,d.dof/2)];
    
end
polyshapes = getPhenotypeFFD(initSamples,d.base);
[~,booleanMap] = getPhenotypeBoolean(polyshapes,d.resolution);

set(0,'DefaultFigureWindowStyle','docked')

for i=1:numel(scale)
    figure(i);
    imagesc(booleanMap{i}); axis equal;
end

axis equal;

%%
tic
[fitness,uMax,uMean,uMin,uStd,data_uMean,data_E] = getFitness(initSamples,d)
toc

%%
set(0,'DefaultFigureWindowStyle','default')

figure(1);
subplot(3,1,1);
plot(data_E,'x-'); 
ylabel('Enstrophy');grid on;
xlabel('Variation \rightarrow');
ax = gca;
ax.XTick = 1:numel(data_E);
ax.XTickLabel = scale(1:end);

subplot(3,1,2);
plot(data_uMean,'x-'); 
ylabel('uMean');grid on;
xlabel('Variation \rightarrow');
ax = gca;
ax.XTick = 1:numel(data_E);
ax.XTickLabel = scale(1:end);

subplot(3,1,3);hold off;
scatter(data_E,data_uMean);
hold on;
X = [ones(length(data_E'),1) data_E'];
b = X\data_uMean'
%b1 = data_E'\data_uMean';
yCalc1 = X*b;
plot(data_E',yCalc1);
grid on;
xlabel('Enstrophy');
ylabel('uMean');

figure(2);
for i=1:numel(booleanMap)
    subplot(4,5,i);
    imagesc(booleanMap{i});
    title(['Variation: ' num2str(scale(i))]);
    axis equal;
end


%%
opts = cmaes;
opts.Restarts = 0;
opts.LBounds = d.ranges(:,1);
opts.UBounds = d.ranges(:,2);
opts.StopFitness = 1e-4;
opts.PopSize = 10;
opts.SaveVariables = 'on';
opts.TolFun = 1e-10;
opts.LogTime = 0;
opts.DispModulo = Inf;
opts.MaxFunEvals = 200;

[XMIN, FMIN, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = cmaes('getFitness', [0.5*ones(1,DOF/2),zeros(1,DOF/2)], [], opts, d)

%%

[~,booleanMap] = getPhenotypeFFD(XMIN',d.base,64);
imagesc(booleanMap{1})
[fitness,uMax,uMean,uMin,uStd] = getFitness(XMIN',d)
%%
[~,booleanMap] = getPhenotypeFFD(BESTEVER.x',d.base,64);
imagesc(booleanMap{1})
%%[fitness,uMax,uMean,uMin,Ustd] = getFitness(BESTEVER.x',d)

%%
his = outcmaesxrecentbest(:,6:end-1)

for i=1:size(his,1)
    his(i,:)
    [~,booleanMap] = getPhenotypeFFD(his(i,:)',d.base,64);
    imagesc(booleanMap{1})
    title(i)
    drawnow
    pause(0.1)
end

%%
[fitness,uMax,uMean,uMin,uStd] = getFitness([0.01*ones(1,DOF/2),zeros(1,DOF/2) ]',d)




