function [value,timings] = mirror_PreciseEvaluate(nextObservations, d)
%mirror_Preci/home/alex/Documents/phd/CFD/snappyhexmesh/mirror/constant/system/controlDictseEvaluate - Send mirror shapes in parallel to OpenFOAM func
%
% Syntax:  [observation, value] = af_InitialSamples(p)
%
% Inputs:
%    nextObservations - [NX1] of parameter vectors
%    d                - domain description struct
%     .openFoamFolder 
%
% Outputs:
%    value(:,1)  - [nObservations X 1] drag force
%
% Other m-files required: mirror_openFoamResult

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: adam.gaier@h-brs.de
% Jun 2017; Last revision: 02-Aug-2017

%------------- BEGIN CODE --------------
folderBaseName = d.openFoamFolder;

% Divide individuals to be evaluated by number of cases
nObs    = size(nextObservations,1);
nCases  = d.nCases;
nRounds = ceil(nObs/d.nCases);
caseStart = d.caseStart;
tic
value = nan(nObs,1);
timings = nan(nObs,1);
for iRound=0:nRounds-1
    PEValue = nan(nCases,1);
    % Evaluate as many samples as you have cases in a batch
    parfor iCase = 1:nCases
    %for iCase = 1:nCases
        disp('NOT running parallel!');
        obsIndx = iRound*nCases+iCase;          
        if obsIndx <= nObs
            openFoamFolder = [folderBaseName 'case' int2str(iCase+caseStart-1) '/']
            PEValue(iCase) = mirror_OpenFoamResult(...
               d.express_CFD(nextObservations(obsIndx,:)),...
               [openFoamFolder 'constant/triSurface/part_07_Mirror.stl'],...
               openFoamFolder, ...
               d.maxDragForce, ...
               d.minDragForce);
        end
    end  
    timings(iRound+1) = toc;
    disp(['Round ' int2str(iRound) ' -- Time so far ' seconds2human(timings(iRound+1))])
    % Assign results of batch 
    obsIndices = 1+iRound*nCases:nCases+iRound*nCases;
    filledIndices = obsIndices(obsIndices<=nObs);
    value(filledIndices) = PEValue(1:length(filledIndices));
end

%------------- END OF CODE --------------