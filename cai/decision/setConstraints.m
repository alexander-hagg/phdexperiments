function constraints = setConstraints(classification,selection,type,method)
%SETCONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

constraints.type                = type;
constraints.selectedClasses     = selection;
constraints.model               = trainConstraintModel(classification.X, classification.simX);

if strcmp(method,'class')
    constraints.classLabels         = classification.labels;
elseif strcmp(method,'individual')
    constraints.classLabels = 1:size(classification.X,1);
end
constraints.threshold           = 0;
constraints.members             = classification.X(ismember(constraints.classLabels, constraints.selectedClasses),:);

end

