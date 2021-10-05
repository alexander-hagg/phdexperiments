function normfeatures = predictFeatures(phenotypes,model)
%GETAEFEATURES Summary of this function goes here
%   Detailed explanation goes here
[~,phenotypes] = getPhenotypeBoolean(phenotypes, model.encoderLG.Layers(1).InputSize(1));

features = getPrediction(phenotypes,model);
if size(features,2) > 2
    features = features*model.pcaCoeff;
end

normfeatures = mapminmax('apply',features',model.normalization);
normfeatures = double(normfeatures');


end

