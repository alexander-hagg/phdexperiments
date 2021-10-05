% see data/basicShape/shapes.m !

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

%% Run cases
clear initSamples;

for i=1:numel(scale)
    sc = scale(i);
    % Small to large circles
    initSamples(i,:) = [sc.*(ones(1,d.dof/2)), zeros(1,d.dof/2)];
    experimentName = 'CirclesSizeVar'
    % Star to circle
    %initSamples(i,:) = [reshape([ones(1,d.dof/4);sc.*ones(1,d.dof/4)],1,d.dof/2), zeros(1,d.dof/2)];
    %experimentName = 'Pointiness'
    
end

% [fitness,polyshapes,booleanMap,movMeanE,movMeanMaxU,stdEnstr,stdU] = getFitness(initSamples,d)

%%
clear enstrophy  maxvelocity enstrophyRunningMean maxvelocityRunningMean
winSize = 1000;
for i=1:length(images)
    cd(int2str(i));
    e = csvread('AllEnstrophies'); ut = csvread('u0'); u = ut(:,1);
    E = movmean(e,min(size(e,1),winSize)); maxU = movmean(u,min(size(e,1),winSize));
    
    enstrophy(i,:) = e;
    enstrophyRunningMean(i,:) = E;
    maxvelocity(i,:) = u;
    maxvelocityRunningMean(i,:) = maxU;
    
    unnormalizedFlowFeatures(i,:) = [E(end),maxU(end)];
    fitness(i) = 1./(1+maxU(end));
    cd('..')
end

%% Check out avg of movmean
    
clear fig;close all;
for i=1:2*length(scale)
    
    fig(i) = figure(i);
    hold off;
    plot(enstrophy(i,:),'k:','LineWidth',1);
    hold on;
    plot(enstrophyRunningMean(i,:),'k-','LineWidth',1);
    plot([1 length(enstrophy(i,:))],[enstrophyRunningMean(i,end) enstrophyRunningMean(i,end)],'b-','LineWidth',1);
    
    if i==1
    legend('Data',['Moving Mean \deltat=' int2str(length(enstrophy(i,:)))],...
        'Final MovMean','Location','NorthEast');
    end
    xlabel('Simulation Steps');
    ylabel('Enstrophy');
    grid on;
    %grid minor;
    ax = gca;
    ax.XTick = 0:500:2000;
    ax.XTickLabel = ax.XTick * 50 ;
    ax.YAxis.Limits = [min(enstrophy(:)) max(enstrophy(:))];
end

%%
save_figures(fig, '.', ['Convergence'], 12, [6 4]);

%%


    