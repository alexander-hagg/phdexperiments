function alignment = evalAlignment(trajectory,varargin)
%EVALALIGNMENT Calculate (single) robot alignment to trajectory. 

dStep = trajectory(2:end,1:2)-trajectory(1:end-1,1:2);
angleTraj = atan2(dStep(:,2),dStep(:,1));
angleRob = trajectory(1:end-1,3);
angleDiff = abs(angleTraj-angleRob);
alignment = mean(angleDiff);



end

