function run_mapelites(numHidden,constraints)
%run_mapelites
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018
%
%------------- BEGIN CODE --------------  

%% Configure experiment
disp("Running MAP-elites on multimodal maze on cluster");
addpath(genpath('.'));
ncores                      = str2num(getenv('NCORES'))-2;
if isempty(ncores); ncores  = 4; end; if isempty(gcp); parpool(ncores); end

d = domain_Maze(numHidden);
p = defaultParamSet(ncores);
p.nGens                     = 2^14;

% Path length
metricFitness               = @(trajectories) sum(sqrt(sum(((trajectories(:,2:end,1:2))-(trajectories(:,1:end-1,1:2))).^2,2)),3);
metricPenalty               = @(x) 1./(x+1e-14); % Fitness is lowest distance, so invert...
evalFcn                     = @(samples) eval_maze(samples,d.numHidden,d.maze,d.timesteps,ncores);
d.fitfun                    = @(X) objective(X, evalFcn, metricFitness, metricPenalty, 'none');

disp(['Running MAP-Elites with ' int2str(d.numHidden) ' hidden nodes; DOF: ' int2str(d.dof)]);
%% Initialize samples
disp('Initialization started');
sobSequence         = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');
sobPoint            = 1;
initSamples         = (2*sobSequence(sobPoint:(sobPoint+p.nChildren)-1,:))-1;
[fitness,trajectories] = d.fitfun(initSamples);
disp('Initialization done');

%% Run MAP-Elites
disp('MAP-Elites started');
    obsMap = createMap(d.featureRes, d.dof);
    [replaced, replacement] = nicheCompete(initSamples, fitness, trajectories, obsMap, d);
    obsMap = updateMap(replaced,replacement,obsMap,fitness,initSamples,[],[]);
    [acqMap, percImproved, percValid, h, allMaps] = mapElites(d.fitfun,obsMap,p,d); 
disp('MAP-Elites done');

save(['out_' int2str(d.numHidden) '.mat'], '-v7.3'); disp('SAVED');

disp('DONE');
end % Function end
