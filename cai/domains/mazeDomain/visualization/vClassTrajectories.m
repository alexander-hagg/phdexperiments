function figHandles = vClassTrajectories(trajectories,values,labels,d,selection,varargin)
%VIEWCLASSTRAJECTORIES Summary of this function goes here
%   labels that are equal to 'nan' are ignored

[uniquelabels] = unique(labels(~isnan(labels))); % remove deselected classes
numClasses = length(uniquelabels);
showVal = 1; if nargin>3; showVal=varargin{1};end

for class=1:numClasses
    %subplot(floor(sqrt(numClasses+1)),ceil(sqrt(numClasses)),class); hold off;
    figHandles(class) = figure;
    classMembers = (labels==uniquelabels(class));
    memberTrajs = trajectories(classMembers,:,:);
    
    vTrajectories(memberTrajs, values(classMembers), d.maze, true, showVal);   
    
    xticks(0:400/d.featureRes(1):400);
    yticks(0:400/d.featureRes(2):400);
    
    addStr = [];
    if showVal == 2 || showVal == 3
        addStr = [int2str(sum(values(classMembers)==1)) '/' int2str(sum(values(classMembers)==2)) '/' int2str(sum(values(classMembers)==3))];        
    end
     if ismember(uniquelabels(class),selection)
         title(['V ' addStr ' Class #' int2str(class)]);
     else
         title(['X ' addStr ' Class #' int2str(class)]);
     end
    hold on;
end

end

