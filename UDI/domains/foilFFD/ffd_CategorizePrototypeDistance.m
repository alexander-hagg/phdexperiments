function [feature] = ffd_CategorizePrototypeDistance(samples, d)
%af_Categorize - Returns feature values between 0 and 1 for each dimension
%
% Syntax:  [feature] = af_Categorize(samples, d)
%
% Inputs:
%   samples  - [N X genomeLength] - uncategorized solutions
%   d        -  Domain description struct
%   .featureMin [1X1]   - minimum feature value
%   .featureMax [1X1]   - maximum feature value
%
% Outputs:
%    feature - [MXN] - Feature values for each individual
%

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 02-Aug-2017

%------------- BEGIN CODE --------------

%% Non vectorized
% for iDeform = 1:size(samples,1)
%    foil = d.express(samples(iDeform,:));
%    [zup(iDeform),xIndx(iDeform)] = max(foil(2,:));
% end
% xup = foil(1,xIndx);

%% Vectorized
model = d.features.constraintModel;
[m_1, s2_1] = gp(model.x.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
    model.trainInput, model.trainOutput(:,1), samples);
[m_2, s2_2] = gp(model.y.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
    model.trainInput, model.trainOutput(:,2), samples);
simspaceCoords = [m_1,m_2];
simspaceCoordConfs = [s2_1,s2_2];

distances = pdist2(simspaceCoords,d.features.prototypes);
feature(:,1) = distances(:,1)./(distances(:,1)+distances(:,2));
feature(:,2) = distances(:,3)./(distances(:,3)+distances(:,4));

%%
feature = (feature-d.featureMin)./(d.featureMax-d.featureMin);
feature(feature>1) = 1; feature(feature<0) = 0;

%------------- END OF CODE --------------
