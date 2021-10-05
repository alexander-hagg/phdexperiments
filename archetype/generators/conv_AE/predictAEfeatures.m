function features = predictAEfeatures(phenotypes,model)
%GETAEFEATURES Summary of this function goes here
%   Detailed explanation goes here
[~,phenotypes] = getPhenotypeBoolean(phenotypes);
features = getPrediction(phenotypes,model);
features = mapminmax('apply',features',model.latentNorm);
features = double(features');

end

