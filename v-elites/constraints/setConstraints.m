function constraints = setConstraints(classification,selection,type,method)
%setConstraints - train constraint model
%
% Syntax:  constraints = setConstraints(classification,selection,type,method)
%
% Inputs:
%    classification
%    selection
%    type
%    method
%
% Outputs:
%    constraints
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Aug 2019; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

constraints.type                = type;
constraints.selectedClasses     = selection;
constraints.model               = trainConstraintModel(classification.X, classification.simX);

if strcmp(method,'class')
    constraints.classLabels     = classification.labels;
elseif strcmp(method,'individual')
    constraints.classLabels     = 1:size(classification.X,1);
end
constraints.threshold           = 0;
constraints.members             = classification.X(ismember(constraints.classLabels, constraints.selectedClasses),:);

end

