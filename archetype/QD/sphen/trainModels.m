function [fitModel,featModels] = trainModels(observation,trueFitness,features,p)
%TRAINMODELS Summary of this function goes here
%   Detailed explanation goes here

fitModel = trainGP(observation,trueFitness,p);
for i=1:size(features,2)
    featModels{i} = trainGP(observation,features(:,i),p);
end

end

