% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% dec 2017; Last revision: 14-Dec-2017

%------------- BEGIN CODE --------------
if nargin < 1
    encoding ='ffd'; nCases=5; startCase=1;
end

% Ensure Randomness of Randomness
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));

% Create Temp Directory for Multithreading
tmpdir = getenv('TMPDIR');
if isempty(tmpdir);tmpdir='/tmp';end
myCluster.JobStorageLocation  = tmpdir;
myCluster.HasSharedFilesystem = true;

%% Add all files to path
addpath(genpath('~/prodigi/'));
cd ~/prodigi;

% Domain
domainname = 'MIRROR'; % domainname = 'FOILFFD';
systemInit;
hpc = true;
if hpc; runFolder = '/home/ahagg2s/prodigi'; else; runFolder = '/home/alex/prodigi'; end
d = mirror_Domain('hpc', hpc, 'lowres', false, 'mirrorOnly', true, 'nCases', nCases, 'runFolder', runFolder);


% Get base and random shape
load('Lea.mat');
observations = mutation;

%% Run OpenFoam
runTime = tic;
[values] = feval(d.preciseEvaluate, observations , d); 
disp(['Runtime: ' seconds2human(toc(runTime))]);

% Save results
save(['foamTest.mat'],'observations', 'values','d');

% set(0,'DefaultFigureWindowStyle','docked')
% 
% horangle = 90;
% angle = 90;
% 
% mutation(1,:) = 0.5*ones(1,51);
% mutation(2,:) = output{1}.prototypes(4,:);
% mutation(3,:) = output{1}.prototypes(11,:);
% 
% for i=1:3
%     figure(i)
%     [FV, ~, ffdP] = mirror_ffd_Express(mutation(i,:), 'mirrorBase.stl');
%     hm = mirrorVisPaper(FV,ffdP, d, horangle, angle, true, false);
%     view(horangle,angle);
%     axis([-120 120 -170 170 600 800]);
% end


