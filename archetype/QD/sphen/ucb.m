function [acquisitionValue, features, predFitness, uncertainty] = ucb(samples, models, d, p)
%acquire - Infill criteria based on uncertainty and fitness
%
% Syntax:  fitness = mirror_AcquisitionFunc(samples, modelPrediction, varCoef)
%
% Inputs:
%   drag -    [2XN]    - drag coefficient mean and variance
%   varCoef  [1X1]    - uncertainty weighting for UCB
%
% Outputs:
%    fitness   - [1XN] - Fitness value (higher is better)
%

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Mar 2020; Last revision: 31-Mar-2020

%------------- BEGIN CODE --------------
fitnessPrediction = predictGP(models(1),samples);

predFitness = fitnessPrediction(:,1);
uncertainty(:,1) = fitnessPrediction(:,2);

for f=2:length(models)
    featurePrediction = predictGP(models(f),samples);
    
    features(:,f-1) = featurePrediction(:,1);
    uncertainty(:,f) = featurePrediction(:,2);
end

if strcmp(p.infill.acqFcn,'FeatureUnCertainty')
    acquisitionValue = (predFitness + ( uncertainty(:,1).*(1+sum(uncertainty(:,2:3)')') )*d.varCoef); 
elseif strcmp(p.infill.acqFcn,'FeatureCertainty')
    acquisitionValue = (predFitness + ( uncertainty(:,1).*(1-sum(uncertainty(:,2:3)')') )*d.varCoef); 
elseif strcmp(p.infill.acqFcn,'PureUCB')
    acquisitionValue = (predFitness + uncertainty(:,1)*d.varCoef); 
end

end

%------------- END OF CODE --------------