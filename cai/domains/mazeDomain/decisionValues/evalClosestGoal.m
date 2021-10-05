function closestGoal = evalClosestGoal(trajectory,goals)
%EVALCLOSESTGOAL Calculate closest goal to end position of a single
%trajectory
%
% TODO: vectorize

[~,closestGoal] = min(sqrt(sum(((goals-trajectory(end,1:2)).^2)')));
end

