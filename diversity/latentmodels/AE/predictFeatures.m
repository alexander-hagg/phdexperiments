function features = predictFeatures(phenotypes,model)
%GETAEFEATURES Summary of this function goes here
%   Detailed explanation goes here
[~,phenotypes] = getPhenotypeBoolean(phenotypes, model.encoderLG.Layers(1).InputSize(1));
features = getPrediction(phenotypes,model);
features = mapminmax('apply',features',model.latentNorm);
features = double(features');

end

