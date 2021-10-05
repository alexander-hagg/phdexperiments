function [fitness,features,polyshapes,booleanMap,unnormalizedFeatures,rawData,allfeatures] = fitnessLettuce(genomes,d, varargin)
%FITNESSLETTUCE Summary of this function goes here
%   [fitness,features,polyshapes,booleanMap] = fitnessLettuce(genomes,d,varargin)

debug = false; if nargin > 2; debug = varargin{1}; end
if debug
    fitness = rand(size(genomes,1),1);
    features = rand(size(genomes,1),2);
    polyshapes = nan(size(genomes,1),1);
    booleanMap = nan(size(genomes,1),1);
    unnormalizedFeatures = features;
    rawData = [];
    return;
end

% Decode to phenotypes
winSize = 1000;
polyshapes = getPhenotypeFFD(genomes,d.base);
[~,booleanMap] = getPhenotypeBoolean(polyshapes,d.resolution);

% Divide individuals to be evaluated by number of cases
nObs    = length(booleanMap);
nCases  = d.nGPUs;
nRounds = ceil(nObs/nCases);

disp(['Number of GPUs: ' int2str(nCases) ' - Number of rounds: ' int2str(nRounds)]);
%%
%tic
for iRound=0:nRounds-1
    disp(['Round ' int2str(iRound+1) '/' int2str(nRounds)]);
    % Evaluate as many samples as you have cases in a batch
    pause(5)
    obsIndices = 1+iRound*nCases:1+(nCases-1)+iRound*nCases;
    parfor iCase = 1:nCases
        obsIndx = iRound*nCases+iCase;
        if obsIndx <= nObs
            workdir = [d.workdir int2str(iCase)];
            system(['rm -r ' workdir '/*']);
            mkdir(workdir);
            system(['cp ' d.projectdir 'domain/footprints/exampleCallLettuce/building.py ' workdir '/']);
            system(['cp ' d.projectdir 'domain/footprints/exampleCallLettuce/buildingDMD.py ' workdir '/']);
            system(['cp ' d.projectdir 'domain/footprints/exampleCallLettuce/DMD_multiple.py ' workdir '/']);
            cd(workdir);
            imwrite(booleanMap{obsIndices(iCase)},'building.png');
            system('touch start');
            system(['cd ' d.workdir]);
        end
    end
    
    for iCase = 1:nCases
        obsIndx = iRound*nCases+iCase;
        if obsIndx <= nObs
            workdir = [d.workdir int2str(iCase)];
            cd(workdir);
            while true
                if exist([workdir '/done'], 'file')
                    e = csvread('AllEnstrophies'); ut = csvread('u0'); u = ut(:,1);
                    E = movmean(e,min(size(e,1),winSize)); maxU = movmean(u,min(size(e,1),winSize));
                    
                    rawData(obsIndx).enstrophy = e;
                    rawData(obsIndx).enstrophyRunningMean = E;
                    rawData(obsIndx).maxvelocity = u;
                    rawData(obsIndx).maxvelocityRunningMean = maxU;
                    
                    unnormalizedFlowFeatures(obsIndx,:) = [E(end),maxU(end)];
                    
                    system('touch stop');
                    system(['cd ' d.workdir]);
                    break;
                else
                    disp('Waiting for case to finish');
                    pause(10);
                end
            end
        end
    end
end
% Normalize features
flowFeatures(:,1) = (unnormalizedFlowFeatures(:,1)-d.featureMin(3))./(d.featureMax(3)-d.featureMin(3));
flowFeatures(:,2) = (unnormalizedFlowFeatures(:,2)-d.featureMin(4))./(d.featureMax(4)-d.featureMin(4));
flowFeatures(flowFeatures>1) = 1; flowFeatures(flowFeatures<0) = 0;

% Shape features
[shapeFeatures,unnormalizedShapeFeatures] = categorize(polyshapes,d);
unnormalizedFeatures = [unnormalizedShapeFeatures,unnormalizedFlowFeatures];

allfeatures = [shapeFeatures,flowFeatures];
features = allfeatures(:,d.selectedFeatures);
fitness = (1./(1+flowFeatures(:,2)))*2-1 ;

cd(d.projectdir);
end

