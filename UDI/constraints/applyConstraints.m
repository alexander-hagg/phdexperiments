function [isValid,classDistances,selectedClasses] = applyConstraints(examples, constraints)
%APPLYCONSTRAINTS Summary of this function goes here
% Detailed explanation goes here
%
% figure(7);hold off
% mmm = ismember(classLabels,selectedClasses);
% scatter(model.trainOutput(mmm,1),model.trainOutput(mmm,2),16,[0 1 0],'filled');hold on
% scatter(model.trainOutput(~mmm,1),model.trainOutput(~mmm,2),16,[1 0 0],'filled');
% scatter(outputs(:,1),outputs(:,2),32,[0 0 0],'filled')
%
% figure(8);hold off;
% tt = model.trainInput(1,:);
% [ott cott] = predictSimspace(tt,model)
% mut = 5*rand(size(tt))-2.5
% [o co] = predictSimspace(tt+mut,model)
% scatter(model.trainOutput(:,1),model.trainOutput(:,2),32,[0.7 0.7 0.7],'filled');hold on;
% scatter(ott(:,1),ott(:,2),32,[0 0 0],'filled');hold on;
% scatter(o(:,1),o(:,2),32,[1 0 0],'filled');
% axis([-10 10 -10 10]);
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 15-Jan-2019

%------------- BEGIN CODE --------------
model = constraints.model; sigma = constraints.threshold; classLabels = constraints.classLabels; selectedClasses = constraints.selectedClasses;
[outputs,confOutputs] = predictSimspace(examples,model);

% Get distances to all class members
for i=1:max(classLabels)
    memberIDs = (classLabels==i);
    
    if constraints.dimreduction
        members = model.trainOutput(memberIDs,:);
        if size(members,1)==1
            classDistances(i,:) = pdist2(outputs,members);
        else
            classDistances(i,:) = min(pdist2(outputs,members)');
        end
    else
        members = model.trainInput(memberIDs,:);
        if size(members,1)==1
            classDistances(i,:) = pdist2(examples,members);
        else
            classDistances(i,:) = min(pdist2(examples,members)');
        end
    end
end

classVariance = sqrt(max(confOutputs'));
classMembership = bsxfun(@le,classDistances-min(classDistances),sigma*classVariance);
isValid = classMembership(selectedClasses,:);
end

