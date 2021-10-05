system(['scp ahagg2s@wr0.wr.inf.h-brs.de:/home/ahagg2s/archetype/trueFitness.csv . && scp ahagg2s@wr0.wr.inf.h-brs.de:/home/ahagg2s/archetype/allMapsSPHEN.mat .']);

load('allMapsSPHEN.mat');
trueFitness = csvread('trueFitness.csv');
nMaps = length(allMaps);


%%
samples = featModels{2}.trainInput;
features = [featModels{1}.trainOutput,featModels{2}.trainOutput];
numReplicates = 500
cvIndices = crossvalind('Kfold',size(samples,1),numReplicates);

%%

for rep=1:numReplicates
    disp([int2str(rep) '/' int2str(numReplicates)]);
    pModel = paramsGP(d.dof);    
    featModel1 = trainGP(samples(cvIndices~=rep,:),features(cvIndices~=rep,1),pModel);
    
    prediction1 = predictGP(featModel1, samples(cvIndices==rep,:));
    trueFeat1 = features(cvIndices==rep,1);    
    mape(1,rep) = mean(abs(prediction1(:,1)-trueFeat1)./mean(trueFeat1))*100;

    featModel2 = trainGP(samples(cvIndices~=rep,:),features(cvIndices~=rep,2),pModel);
    prediction2 = predictGP(featModel2, samples(cvIndices==rep,:));
    trueFeat2 = features(cvIndices==rep,2);    
    mape(2,rep) = mean(abs(prediction2(:,1)-trueFeat2)./mean(trueFeat2))*100;

end

%%
disp(mape)
close all;

figure(1);hold off;
plot(prediction2(:,1))
hold on;
plot(trueFeat2)
%plot(mape)


