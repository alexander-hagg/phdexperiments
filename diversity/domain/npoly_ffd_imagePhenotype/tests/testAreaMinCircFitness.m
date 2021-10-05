function testAreaMinCircFitness(matchedParams,d)
%TESTAREAMINCIRCFITNESS Summary of this function goes here
%   Detailed explanation goes here
%% fitness test

load('matchedParams.mat');
[fitness,phenotypes]  = d.fitfun(matchedParams)
figHandle = figure(100);
for i=1:length(phenotypes)
    subplot(2,2,i);
    showPhenotype(matchedParams(i,:),d,figHandle);axis([-1 1 -1 1]);
    title(['Fitness: ' num2str(fitness(i))]);
end

end

