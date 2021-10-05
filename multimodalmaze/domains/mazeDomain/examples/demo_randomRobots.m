%demo_randomRobots
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018
%
%------------- BEGIN CODE --------------

clear
maze = 'mediumRound'; %'exampleRound'
ncores = 1;
numHidden = 2;
numNetworks = 5;

% Configure experiment
d = domain_Maze(numHidden,ncores);

sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');
initSamples = 1*((2*sobSequence(1:(1+numNetworks)-1,:))-1);

timesteps = 5000;

d.debug = false;
[trajectories,values] = eval_maze(initSamples,d.numHidden,d.maze,timesteps,d.ncores,d.goalLocations,d.useRNN,d.debug,false,d.useGoalQuadrantSensors);
%
plottingMap = true;
figure(1);hold off;
for i=1:size(trajectories,1)
    trajectory = squeeze(trajectories(i,:,:));
    plotSim(trajectory,d,plottingMap);hold on;
    plottingMap = false;
end
%initSamples = initSamples + 0.01*randn(size(initSamples));
%%
figure(2);hold off;
[~,dist] = evalTrajectoryClosestToGoal(trajectory,d.goalLocations);
plot(dist(1,:)','Color',[1 0 0],'LineWidth',4);hold on;
plot(dist(2,:)','Color',[0 1 0],'LineWidth',4);
plot(dist(3,:)','Color',[0 0 1],'LineWidth',4);
plot(25*ones(1,size(dist,2)),':','Color',[0 0 0],'LineWidth',4);
legend('Goal 1','Goal 2','Goal 3');
title('Distance to goal');axis([0 size(dist,2) 0 300]);


