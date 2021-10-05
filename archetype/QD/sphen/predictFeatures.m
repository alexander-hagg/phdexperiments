function [features,confidence] = predictFeatures(x,models)
%PREDICTFEATURES Predict features using GP models
%
% Syntax:  features = predictFeatures(x,models)
%
% Inputs:
%    x - locations (genomes)
%    models - trained GP models, one for each feature
%
% Outputs:
%    features - predicted features
%
%

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Apr 2020; Last revision: 03-Apr-2020

for m=1:length(models)
     feat = predictGP(models{m},x);
     features(:,m) = feat(:,1);
     confidence(:,m) = feat(:,2);
end
end

