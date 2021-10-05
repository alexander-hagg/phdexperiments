function [outputs,confOutputs] = predictSimspace(locations,model,confidenceAsDim)
%predictSimspace - predict locations in similarity space
%
% Syntax:  [outputs,confOutputs] = predictSimspace(locations,model)
%
% Inputs:
%    locations        - [NxD] - N samples with dimensionality D.
%    model            - similarity space model
%
% Outputs:
%    outputs           - [Nx2] - predicted locations
%    confOutputs       - [Nx2] - model confidences of locations
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

[m_1, s2_1] = gp(model.x.hyp, @infExact, model.x.meanfunc, model.x.covfunc, model.x.likfunc, ...
    model.trainInput, model.trainOutput(:,1), locations(~any(isnan(locations')),:));
[m_2, s2_2] = gp(model.y.hyp, @infExact, model.y.meanfunc, model.y.covfunc, model.y.likfunc, ...
    model.trainInput, model.trainOutput(:,2), locations(~any(isnan(locations')),:));
outputs = nan(size(locations,1),2);
outputs(~any(isnan(locations')),:) = [m_1,m_2];

confOutputs = nan(size(locations,1),2);
confOutputs(~any(isnan(locations')),:) = [s2_1,s2_2];

outputs = outputs';
confOutputs = confOutputs';
if confidenceAsDim
    outputs = [outputs;mean(confOutputs)];
end
end

%------------- END CODE --------------
