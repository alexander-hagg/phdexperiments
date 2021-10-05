function [closenessID,dist] = evalTrajectoryClosestToGoal(trajectory,goals)
%EVALCLOSESTGOAL Calculate closest goal to end position of a single
%trajectory
%
% TODO: vectorize

for i=1:size(goals,1)
    dist(i,:) = sqrt(sum((goals(i,:)-trajectory(:,1:2)).^2,2));    
end

[~,id] = min(dist(:));
[row, col] = ind2sub(size(dist),id);
closenessID = row;
end

