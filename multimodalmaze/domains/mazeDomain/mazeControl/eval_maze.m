function [trajectories,values] = eval_maze(weights,d,varargin)
%eval_maze
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018
%
%------------- BEGIN CODE -------------- 
visSim = false; if nargin > 2; visSim = varargin{1};end

numEvals = size(weights,1);
runs = ceil(numEvals./d.ncores);
allIDs = 1:numEvals;

if ~d.useRNN
    simCmd = 'domains/mazeDomain/simulator/build/./fastsim-custom';
else
    simCmd = 'domains/mazeDomain/simulator/build/./fastsim-custom-rnn';
end

trajectories = zeros(numEvals,d.timesteps,3);
for p = 1:runs
    ID = ((p-1)*d.ncores+1):p*d.ncores;
    ID = ID(ismember(ID,allIDs));
    for iter = 1:length(ID)
        thisID = ID(iter); weightsStr = [];
        for i=1:size(weights(thisID,:),2)
            weightsStr = [weightsStr ' ' num2str(weights(thisID,i))];
        end
        trajFilename = [d.tmpdir  '/traj.dat' int2str(iter)];
        mkdir(d.tmpdir);
        if visSim
            callEval = [simCmd ' domains/mazeDomain/simulator/worlds/' d.mazeCfgFile 'Vis.xml ' d.tmpdir  ' ' int2str(d.timesteps) ' ' int2str(iter) ' ' int2str(d.useGoalQuadrantSensors) ' ' int2str(d.numHidden) ' ' int2str(d.debug) weightsStr];
        else
            callEval = [simCmd ' domains/mazeDomain/simulator/worlds/' d.mazeCfgFile '.xml ' d.tmpdir  ' ' int2str(d.timesteps) ' ' int2str(iter) ' ' int2str(d.useGoalQuadrantSensors) ' ' int2str(d.numHidden) ' ' int2str(d.debug) weightsStr];
        end
        system(callEval);

        if exist(trajFilename, 'file')
            data = dlmread(trajFilename,' ');
        
            trajectories(thisID,1:size(data,1),1:3) = data(:,1:3);
            realTraj = squeeze(trajectories(thisID,:,1:3));
            [~,ids] = unique(realTraj,'rows'); ids = sort(ids); realTraj = realTraj(ids,:);

            % Calculate robot alignment to trajectory
            values{1}(thisID) = evalAlignment(realTraj);        
            % Calculate first exit robot took on selected ring 
            for i=1:length(d.goalRings)
                values{i+1}(thisID) = evalFirstGoal(realTraj,d.goalRings{i},d.ringDiameter(i),d.center);
            end
        else
            % Just leave trajectory to zeros and set values to 0 as well
            for i=1:4
                values{i}(thisID) = 0;
            end
        end
       
    end
    clear theseTrajectories
end

end

