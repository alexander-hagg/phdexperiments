function [classDistances,simspacePredictions] = applyConstraints(samples, c)
%applyConstraints - apply user selection constraints
%
% Syntax:  [isValid,classDistances,selectedClasses,outputs] = applyConstraints(samples, constraints)
%
% Inputs:
%    samples        - [NxD] - N samples with dimensionality D
%    c              - struct - user constraints
%
% Outputs:
%    classDistances      - distances to closest members of classes
%    simspacePredictions - predicted similarity space locations of samples
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------
simspacePredictions = predictSimspace(samples,c.model);

% Get distances to all class members
for i=1:max(c.classLabels)
    memberIDs = (c.classLabels==i);
    members = c.model.trainOutput(memberIDs,:);
    if size(members,1)==1
        classDistances(i,:) = mypdist2(simspacePredictions,members);
    else
        classDistances(i,:) = min(mypdist2(simspacePredictions,members)');
    end
end

end

