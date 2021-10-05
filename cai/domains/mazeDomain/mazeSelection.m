function [isSelected,selection] = mazeSelection(estimatedLabels, values, critID, critVal, critThresh)
%SELECT Summary of this function goes here
%   Detailed explanation goes here
selection = [];
classlabels = unique(estimatedLabels);
for ii=1:max(classlabels)
    if critID == 1
        disp('Selection procedure for criterion 1, wiggliness of the path, does not work yet!');
        %if critVal<=0.5
        %    percCorrect = sum(values{critID}(estimatedLabels==ii)<critVal)/length(values{critID}(estimatedLabels==ii));
        %elseif critVal>=0.5
        %    percCorrect = sum(values{critID}(estimatedLabels==ii)>critVal)/length(values{critID}(estimatedLabels==ii));
        %end
        %if percCorrect > critThresh
        %    selection(end+1) = ii;
        %end
    elseif critID > 1
        comparedVals = values{critID}(estimatedLabels==ii);
        compareWith = [0 1 2 3];          % Get number of exits in the maze indirectly from configuration
        a = hist(comparedVals,compareWith);         % Count
        if (a(critVal+1)/sum(a(2:end))) >= critThresh
            selection(end+1) = ii;
        end
    else
        disp('No valid selection criterion, or for this selection criterion, no procedure was programmmed yet');
    end
end
isSelected = ismember(estimatedLabels,selection);
end
