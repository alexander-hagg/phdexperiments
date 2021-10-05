function [outputs,confOutputs] = predictSimspace(locations,model)
%PREDICTSIMSPACE Summary of this function goes here
%   Detailed explanation goes here
[m_1, s2_1] = gp(model.x.hyp, @infExact, model.x.meanfunc, model.x.covfunc, model.x.likfunc, ...
    model.trainInput, model.trainOutput(:,1), locations(~any(isnan(locations')),:));
[m_2, s2_2] = gp(model.y.hyp, @infExact, model.y.meanfunc, model.y.covfunc, model.y.likfunc, ...
    model.trainInput, model.trainOutput(:,2), locations(~any(isnan(locations')),:));
outputs = nan(size(locations,1),2);
outputs(~any(isnan(locations')),:) = [m_1,m_2];

confOutputs = nan(size(locations,1),2);
confOutputs(~any(isnan(locations')),:) = [s2_1,s2_2];
end

