clear;clc;
DOF = 16;DOMAIN = 'npoly_ffd_imagePhenotype';
addpath(genpath('/home/alex/poem'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));

d = domain(DOF);
p = poemParamSet;
%d.fitfun = d.fitfunAreaCirc; % Unimodal function
d.fitfun = d.fitfunPointSymm; % Multimodal function

m = getAEConfig('data/workdir');

sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,phenotypes] = d.fitfun(initSamples);

%%
[map,m,stats] = poem(initSamples,fitness,phenotypes,p,d,m);


%%
figure(5);hold off;
plot(median(stats.model.trainingLosses'));
title('Median AE training loss');

figure(6);hold off;
plot(stats.error.mean);hold on;
plot(stats.error.median);
axis([1 length(stats.error.mean) 0 1000]);title('Reconstruction Error Map');legend('Mean','Median');


figure(7);hold off;
plot(stats.fitness.mean);hold on;
plot(stats.fitness.median);
axis([1 length(stats.fitness.mean) 0 1]);title('Fitness Map');legend('Mean','Median');

figure(8);hold off;
plot(stats.fitness.total);hold on;
%axis([1 length(stats.fitness.mean) 0 1]);
title('Total Fitness Map');

figure(9);hold off;
plot(stats.elites.number);hold on;
%axis([1 length(stats.fitness.mean) 0 1]);
title('Total Number of Solutions');

%figure(9);hold off;
%plot(mean(stats.model.trainingLosses));hold on;
%plot(median(stats.model.trainingLosses));hold on;
%plot(std(stats.model.trainingLosses));hold on;
%axis([1 length(mean(stats.model.trainingLosses)) 0 1000]);
%title('Avg. Training Error AE');

%%
placement = 64*stats.maps{end}.features;
figHandle = figure(11); showPhenotype(stats.maps{end}.genes,d,figHandle,placement);
title('Polygons maximizing area/perimeter');

%% Show polygon maps evolving

for i=[1 ceil(length(stats.maps)/2) length(stats.maps)]
    figHandle = figure(i);
    placement = 64*stats.maps{i}.features;
    showPhenotype(stats.maps{i}.genes,d,figHandle,placement);    
    drawnow;
end




%% END OF CODE