function [goalID,dist] = evalFirstGoal(trajectory,goals,diameter,center)
%EVALFIRSTGOAL Return first goal reached within threshold
%   Ignore all trajectory values outside of the maze
%   Radius of maze = 336 pixels

threshold   = 30;
tolerance   = 2;
goals(4,:)  = center;     % Measure distance to origin to exclude trajectories outside of the maze
dist = pdist2(goals,trajectory(:,1:2));


outCircle = 2*(dist(4,:)>((diameter/2)+tolerance));
inCircle = (dist(4,:)<((diameter/2)-tolerance));

outCircle = outCircle + inCircle; % 2 is out, 1 is in, 0 is in tolerance zone

lastPointIncircle = find(fliplr(outCircle)==1); lastPointIncircle = lastPointIncircle(1);
lastPointIncircle = size(trajectory,1)-lastPointIncircle; % flip back!

%figure(1);hold off;plot((dist(4,:)>(diameter/2)),'LineWidth',4); hold on;
%scatter(lastPointIncircle,0,'filled');

outCircle = zeros(1,length(outCircle));
outCircle(1:lastPointIncircle) = 0;
outCircle(lastPointIncircle+1:end) = 1;

%plot(outCircle,'LineWidth',4); hold on;

goalIDs(1,:) = ((dist(1,:)<threshold)') & outCircle';
goalIDs(2,:) = ((dist(2,:)<threshold)') & outCircle';
goalIDs(3,:) = ((dist(3,:)<threshold)') & outCircle';

[firstGoalReached,col] = find(goalIDs); 

if isempty(firstGoalReached)
    goalID = 0;
else
    goalID = firstGoalReached(1);
end
end

