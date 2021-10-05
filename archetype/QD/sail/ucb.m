function [acquisitionValue, features, predFitness, uncertainty] = ucb(samples, modelPrediction, d, p)
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
fitnessPrediction = predictGP(modelPrediction,samples);
predFitness = fitnessPrediction(:,1);
uncertainty = fitnessPrediction(:,2);
% better fitness is higher fitness  
acquisitionValue = (predFitness + (uncertainty*d.varCoef)); 

phenotypes = d.getPhenotype(samples);
features = d.categorize(samples,phenotypes,p,d);

end

%------------- END OF CODE --------------