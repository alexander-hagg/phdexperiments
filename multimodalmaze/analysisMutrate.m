function analysisMutrate(folder)
%
% folder                - Should contain a folder name (char array)
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2019; Last revision: 30-Jan-2019
%
%------------- BEGIN CODE --------------

addpath(genpath('.'));

d = domain_Maze;
angles = [30 90 150 210 270 330];
numA =length(angles);
%mutationRates = [0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.015 0.02 0.05 0.08 0.10 0.15 0.20];

% Planning
mutationRates = [0.0001 0.0005 0.001 0.005 0.01 0.015 0.02 0.05 0.08];

% Controller 
%mutationRates = [0.00001 0.0001 0.001 0.01 0.02 0.05 0.10 0.20];

methods = {'UDHM','SEED','COMB'};
exits = [1 2 3];
load([folder '/' 'selectedMaps.mat']);

fdataRunData = dataRunData;
fdataRunData(cellfun(@isempty,{dataRunData(:,:,:,:).runData})) = [];
fdataRunData = reshape(fdataRunData,length(methods),length(mutationRates),length(angles),[]);

disp('Data loading complete');
%%
clear exitCorrectness exitInCorrectness mCor mInCor mCor25 mCor75 mInCor25 mInCor75 userSelDrift
countMaps = zeros(2,1);
for methodID=1:length(methods)
    disp(['Running for method ' int2str(methodID) '/' int2str(length(methods))]);
    for mutRateID=1:length(mutationRates)
    disp(['Running for mutation rate ' int2str(mutRateID) '/' int2str(length(mutationRates))]);
        
        for angle=1:numA
            disp('switch angle');
            for exitID=[1:3]
                disp('switch exit');
                %dataO = fdata(methodID,mutRateID,angle,exitID);
                dataA = fdataRunData(methodID,mutRateID,angle,exitID);
                if isempty(dataA.runData); continue; end
                exitCorrectness(methodID,mutRateID, (angle-1)*(3) + (exitID)  ) = sum(sum(dataA.runData{end}.acqMap.ring1(:)==dataA.runData{end}.p.selectionValue));
                exitInCorrectness(methodID,mutRateID, (angle-1)*(3) + (exitID)  ) = sum(sum(dataA.runData{end}.acqMap.ring1(:)~=dataA.runData{end}.p.selectionValue)) - sum(isnan(dataA.runData{end}.acqMap.fitness(:)));
                
                if methodID==2
                    X = reshape(dataA.runData{end}.acqMap.genes,900,size(dataA.runData{end}.acqMap.genes,3));
                    X = X(all(~isnan(X')),:);
                    penalty = constraintPenalty(X,dataA.runData{2}.constraints);
                    userSelDrift(methodID,mutRateID, (angle-1)*(3) + (exitID) ) = nanmedian(penalty(:));
                    
                else
                    userSelDrift(methodID,mutRateID, (angle-1)*(3) + (exitID) ) = nanmedian(dataA.runData{end}.acqMap.penalty(:));
                end
            end
        end
    end
end

%%
normFactorCorrectness = 900;

mCor = 100*median(exitCorrectness./normFactorCorrectness,3);
mCor25 = 100*prctile(exitCorrectness./normFactorCorrectness,25,3);
mCor75 = 100*prctile(exitCorrectness./normFactorCorrectness,75,3);

mInCor = 100*median(exitInCorrectness./normFactorCorrectness,3);
mInCor25 = 100*prctile(exitInCorrectness./normFactorCorrectness,25,3);
mInCor75 = 100*prctile(exitInCorrectness./normFactorCorrectness,75,3);

mDrift = median(userSelDrift,3);
mDrift25 = prctile(userSelDrift,25,3);
mDrift75 = prctile(userSelDrift,75,3);




save([folder '/' 'analysis.mat'],'mCor','mCor25','mCor75','mInCor','mInCor25','mInCor75','mDrift','mDrift25','mDrift75', ...
    'mutationRates', 'methods', 'angles', 'exitCorrectness', 'exitInCorrectness', 'userSelDrift', 'normFactorCorrectness');

disp('SAVED and DONE');

end

