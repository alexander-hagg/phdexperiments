clear;clc;
load('initsamples64.mat');
franges = [movMeanE;movMeanMaxU];
shapeFeatures = d.categorize([],polyshapes,[],d)';

load('pointiness.mat');
franges2 = [movMeanE;movMeanMaxU];
shapeFeatures2 = d.categorize([],polyshapes,[],d)';

load('circles.mat');
franges3 = [movMeanE;movMeanMaxU];
shapeFeatures3 = d.categorize([],polyshapes,[],d)';

features = [franges,franges2,franges3;shapeFeatures,shapeFeatures2,shapeFeatures3];
featureNames = {'E','uMax','Area','Perimeter'};
% 64 9 9
%%
set(0,'DefaultFigureWindowStyle','default')

clear figs;
featIDs = 1:4;
c = combvec(featIDs,featIDs);

for i=1:size(c,2)
    if c(1,i) == c(2,i)
        continue;
    end
fig(i) = figure(i);hold off;
scatter(features(c(1,i),1:64),features(c(2,i),1:64),32,'k');
hold on;
scatter(features(c(1,i),65:73),features(c(2,i),65:73),32,'r','filled');
scatter(features(c(1,i),74:82),features(c(2,i),74:82),16,'b','filled');
ax = gca; 
%axis([0.4 2 0.08 0.18]);
axis([min(features(c(1,i),:)) max(features(c(1,i),:)) min(features(c(2,i),:)) max(features(c(2,i),:))]);
legend('Random Initial Set','Pointiness Set', 'Circles Set');
xlabel(featureNames{c(1,i)});
ylabel(featureNames{c(2,i)});

grid on;

end
%%

save_figures(fig, '.', ['lettuceFeatureRanges'], 12, [5 4]);
