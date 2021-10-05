function isValid = applyConstraints(examples, model, sigma, classLabels, selectedClasses)
%APPLYCONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

[m_1, s2_1] = gp(model.x.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
    model.trainInput, model.trainOutput(:,1), examples);
[m_2, s2_2] = gp(model.y.hyp, @infExact, model.p.meanfunc, model.p.covfunc, model.p.likfunc, ...
    model.trainInput, model.trainOutput(:,2), examples);
outputs = [m_1,m_2];
confOutputs = [s2_1,s2_2];

%% Get sets of candidate classes

% Get distances to all class members
for i=1:max(classLabels)
    memberIDs = (classLabels==i);
    members = model.trainOutput(memberIDs,:);
    % Get distance to class for all outputs
    deltaClasses(i,:) = min(pdist2(outputs,members)');
end
delta = deltaClasses - min(deltaClasses); % Closest class at zero
variance = sqrt(max(confOutputs'));
classMembership = bsxfun(@le,delta,sigma*variance);

isValid = classMembership(selectedClasses,:);
end

