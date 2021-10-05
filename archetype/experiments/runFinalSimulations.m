clear;clc;
DOF = 16;
workpath = ['/home/' getenv('USER') '/archetype/'];
addpath(genpath(workpath));
DOMAIN = 'footprints';  rmpath(genpath([workpath 'domain'])); addpath(genpath([workpath 'domain/' DOMAIN]));
QD = 'grid';            rmpath(genpath([workpath 'QD/grid'])); rmpath(genpath([workpath 'QD/voronoi'])); addpath(genpath([workpath 'QD/' QD]));
SAQD = 'sphen';          rmpath(genpath([workpath 'QD/sail'])); rmpath(genpath(['QD/sphen'])); addpath(genpath([workpath 'QD/' SAQD]));

nGPUs = str2num(getenv('numGPUs')); if isempty(nGPUs); nGPUs = 4; end

load('finalResult.mat');
%%
nObs    = size(genes,1);
nCases  = nGPUs;
nRounds = ceil(nObs/nCases);

saveDir = [d.workdir 'results/'];
mkdir(saveDir);

for iRound=0:nRounds-1
    disp(['Round ' int2str(iRound+1) '/' int2str(nRounds)]);
    obsIndices = 1+iRound*nCases:1+(nCases-1)+iRound*nCases;
    genomeRanges = obsIndices(obsIndices <= nObs);    
    fitnessLettuce(genes(genomeRanges,:),d);
    
    for i=1:length(genomeRanges)
        system(['cp -r ' d.workdir int2str(i) ' ' saveDir int2str(genomeRanges(i))]);
    end
end