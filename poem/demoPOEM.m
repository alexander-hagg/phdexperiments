clear;clc;
DOF = 16;DOMAIN = 'npoly_ffd_imagePhenotype';
addpath(genpath('/home/alex/poem'));rmpath(genpath('domain')); addpath(genpath(['domain/' DOMAIN]));

d = domain(DOF);
p = poemParamSet;
m = getAEConfig('data/workdir');

sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,phenotypes] = d.fitfun(initSamples);

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
%%
placement = 64*stats.maps{end}.features;
figHandle = figure(8); showPhenotype(stats.maps{end}.genes,d,figHandle,placement);
title('Autodiscovered Manifold of Point Symmetric Polygons');


%% END OF CODE