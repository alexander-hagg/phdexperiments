function [isSelected,runData] = selectionProcedure(values, classification, p, runData)
%SELECTIONPROCEDURE Selection procedure for maze domain
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2019; Last revision: 03-Jan-2019
%
%------------- BEGIN CODE --------------
disp(['Selecting started: selecting classes based on ' p.selectionMethod]);
thresholds = 1:-0.01:0.01;

for selThresID = 1:length(thresholds)
    selThresh = thresholds(selThresID);
    if strcmp(p.selectionMethod,'class')
        labels = classification.labels;
    elseif strcmp(p.selectionMethod,'individual')
        labels = 1:length(values{2});
    end
    [isSelected,selection] = mazeSelection(labels, values, ...
        p.selectionCriterionID, p.selectionValue, selThresh);
    if (selThresh < p.selectionThreshold) && sum(isSelected) > p.minSelectedSamples
        disp(['Threshold at: ' num2str(selThresh)]);
        runData.selThresh = selThresh;
        break;
    end
end

if sum(isSelected) < p.minSelectedSamples
    disp('Not enough samples selected, skipping selection process for now');
    runData.constraints = 'none';
else
    runData.isSelected = isSelected; runData.selection = selection;
    %% Train Constraints
    constraints = setConstraints(classification, selection, p.constraintType, p.selectionMethod);
    constraints.dimreduction = p.useDimReduction;
    disp(['Using dimreduction: ' mat2str(logical(constraints.dimreduction))]);
    runData.constraints = constraints;
    disp(['Done training constraint model']);
    
end

end
%------------- END CODE --------------