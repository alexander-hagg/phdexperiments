function [trajectories,values] = eval_maze_plan(weights,d,varargin)
%eval_maze
% [trajectories,values] = eval_maze_plan(weights,d,varargin)
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2018; Last revision: 17-Jan-2018
%
%------------- BEGIN CODE --------------
numSamples  = size(weights,1);

% Calculate phenotype
for i=1:numSamples
    phenotype(i,:,:) = phenoPlan(weights(i,:),d);
end
    
trajectories = zeros(numSamples,d.phenotypeLength,3);
for i=1:numSamples
    if any(isnan(phenotype(i,:)))
        trajectories(i,:,:) = nan(d.phenotypeLength,3);
        for ringID=1:length(d.goalRings)
            values{ringID+1}(i) = nan;
        end        
        continue;
    end
    [cx,cy,c] = improfile(d.map,phenotype(i,1,:),phenotype(i,2,:),d.phenotypeLength);
    
    trajectories(i,:,1:2) = [cx,cy]; %[cx,cy,zeros(length(cx),1)];
    realTraj = squeeze(trajectories(i,:,:));
    %[~,ids] = unique(realTraj,'rows'); ids = sort(ids); realTraj = realTraj(ids,:);
    
    % Calculate robot alignment to trajectory
    values{1}(i) = 0; %%evalAlignment(realTraj);
    % Calculate first exit robot took on selected ring
    for ringID=1:length(d.goalRings)
        values{ringID+1}(i) = evalFirstGoal(realTraj,d.goalRings{ringID},d.ringDiameter(ringID),d.center);
    end
end




end

