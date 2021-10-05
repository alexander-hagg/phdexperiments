function model = trainAEfeatures(phenotypes,model)
%TRAINAEFEATURES Summary of this function goes here
%   Detailed explanation goes here
[~,phenotypes] = getPhenotypeBoolean(phenotypes);
model = model.train(phenotypes);
[~,latentNorm] = mapminmax(model.latent,0,1);
model.latentNorm = latentNorm;

end

