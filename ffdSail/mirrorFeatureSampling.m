%PRODIGI vs SAIL comparison -
%
% Other submodules required: sail, gpml
%

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2018; Last revision: 25-Jan-2018

%------------- BEGIN CODE --------------
% Clean up workspace and add relevant files to path
clear;
currentPath = mfilename('fullpath');
addpath(genpath(currentPath(1:end-length(mfilename))));

% Domain
domainname = 'MIRROR';%domainname = 'FOILFFD';
systemInit;
d = mirror_Domain('hpc', false, 'lowres', true); %d = ffd_Domain; %d = velo_Domain; %d = af_Domain;
%d = ffd_Domain('/scratch/ahagg2s/sailCFD/tmp');

nInitialSamples = 10

% Produce initial solutions
sobSequence  = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');
sobPoint = 1;

% Replace with next in Sobol Sequence
nextGenes = sobSequence(sobPoint:(1+nInitialSamples),:);

d.express(nextObservations(obsIndx,:))