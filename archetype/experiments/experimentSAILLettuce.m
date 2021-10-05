%experimentSAILLettuce - Run Lettuce CFD domain
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Apr 2020; Last revision: 11-Apr-2019
%
%------------- BEGIN CODE --------------

clear;clc;
DOF = 16;
workpath = ['/home/' getenv('USER') '/archetype/'];
addpath(genpath(workpath));
DOMAIN = 'footprints';  rmpath(genpath([workpath 'domain'])); addpath(genpath([workpath 'domain/' DOMAIN]));
QD = 'grid';            rmpath(genpath([workpath 'QD/grid'])); rmpath(genpath([workpath 'QD/voronoi'])); addpath(genpath([workpath 'QD/' QD]));
SAQD = 'sail';          rmpath(genpath([workpath 'QD/sail'])); rmpath(genpath(['QD/sphen'])); addpath(genpath([workpath 'QD/' SAQD]));

nGPUs = str2num(getenv('numGPUs')); if isempty(nGPUs); nGPUs = 4; end
d = domain(DOF,'lettuce_SAIL',nGPUs);
p = defaultParamSet;

p.infill = infillParamSet;
p.infill.nTotalSamples = 1024;
p.infill.nAdditionalSamples = 16;

sailP = p; sailP.nGens = 2^12;  sailP.nChildren = 2^5;

d.fitfun = d.fitfunLettuce; % Multimodal function: point symmetry (center); for testing purposes
experimentName = 'LETTUCE_SAIL'
numReplicates = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUGGING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d.fitfunLettuce = @(geno) fitnessLettuce(geno,d,true);
% d.fitfun = d.fitfunLettuce;
% sailP = p; sailP.nGens = 2^7;  sailP.nChildren = 2^5;
% sailP.infill.nTotalSamples = 80
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /DEBUGGING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


initSampleSet = 'initsamples64.mat';
if exist(initSampleSet,'file')
    load('initsamples64.mat');
    initSet.features = initSet.features(:,d.selectedFeatures);
else    
    disp('Creating initial sample set');
    sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen'); sobPoint = 1;
    initSet.samples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+p.numInitSamples)-1,:)+d.ranges(:,1)';
    [initSet.fitness,initSet.features] = d.fitfun(initSet.samples);
    save('initsamples64.mat','initSet');
end
%% ----------------------------------------------------------------------------------

for rep=1 : numReplicates
    %% SAIL
    [mapSAIL{rep},surrogateFitnessSAIL{rep},allMapsSAIL{rep}] = sail(initSet,sailP,d,1);    
    save([experimentName '_replicate_' int2str(rep) '.mat']);
end



