function dataExtractMutrate(folder, numStepsToExtract, recursion)
%%DATAEXTRACTMULTIRING extract partial experimental data for multiring multimodal maze experiment with optional recursion
%
% folder                - Should contain a folder name (char array)
% numStepsToExtract     - Number of time steps per experiment blob
%                           if == 1, extract last step
%                           if == 2, extract first and last step
%                           if == n, extract n steps, including first and
%                           last
% recursion             - recurse over folders
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2019; Last revision: 09-Jan-2019
%
%------------- BEGIN CODE --------------

%%
addpath(genpath('.'));
eval(['delete ' folder '/' 'selectedMaps.mat']);

% Collect folders recursively
if recursion; folders = getFoldersRecursive(folder); end

angles          = [30 90 150 210 270 330];
numA            = length(angles);
% mutationRates   = [0.000001 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.015 0.02 0.05 0.08 0.10 0.15 0.20];
% mutationRates   = [0.00001 0.0001 0.001 0.01 0.02 0.05 0.10 0.20];
mutationRates   = [0.00001 0.0001 0.001 0.005 0.01 0.015 0.02 0.05 0.08];
methods         = {'UDHM','SEED','COMB'};
exits           = [1 2 3];



for fr=1:length(folders)
    % Find 'allMaps' files, which contain each iteration's data
    fsPaths = dir(folders{fr});
    fsPaths([fsPaths.isdir]) = [];
    wanted = [];
    for fi=1:length(fsPaths)
        k = strfind(fsPaths(fi).name,'allMaps');
        wanted(fi) = ~isempty(k);
    end
    fsPaths(~wanted) = [];
    
    % Open files and extract maps and description
    for fi=1:length(fsPaths)
        fname = [folders{fr} '/' fsPaths(fi).name]
        
        %% Deduce experimental parameters
        angle = regexp(fsPaths(fi).name,'\d*_allMaps.mat','Match'); angle = str2num(angle{1}(1:end-12));
        angleID = find(angle==angles);
        
        exitIDs = strfind(folders{fr},'exit');
        if ~isempty(exitIDs)
            exitNR = str2num(folders{fr}(exitIDs(1)+4));
        else
            continue;
        end
        
        methods = {'UDHM','SEED','COMB'};
        isUDHM = strfind(folders{fr},'UDHM');
        if ~isempty(isUDHM), isUDHM = true; else, isUDHM = false;end
        isSEED = strfind(folders{fr},'SEED');
        if ~isempty(isSEED), isSEED = true; else, isSEED = false;end
        
        isCOMB = strfind(folders{fr},'COMB');
        if ~isempty(isCOMB)
            isSEED = true; isUDHM = true;
        else
            isCOMB = false;end
        
        if  (~isSEED && ~isUDHM)
            disp('Uhh.. experiment is neither UDHM or SEED? Let us ignore this one then.. *slowly moves to the next iteration*'); 
            continue;
        else
            methodID = isUDHM*1 + isSEED*2;            
        end
        
        mutRate = regexp(fsPaths(fi).folder,'0\.\d*','Match'); mutRate = str2num(mutRate{1});
        mutRate = find(mutRate==mutationRates);
        %%
        
	
        data(methodID,mutRate,angleID,exitNR).fname = fname;
        data(methodID,mutRate,angleID,exitNR).description.mutRate = mutRate/1000;
        data(methodID,mutRate,angleID,exitNR).description.angle = angle;
        data(methodID,mutRate,angleID,exitNR).description.methodID = methods{methodID};
        % load(fname,'runData_EXTRA');
        load([fsPaths(fi).folder '/' int2str(angle) '.mat'],'runData');
        dataRunData(methodID,mutRate,angleID,exitNR).runData = runData;
        
        % for iter=1:length(runData_EXTRA)
        %    if isempty(runData_EXTRA{iter}); continue; end
        %    data(methodID,mutRate,angle,exitNR).maps(iter,:) = runData_EXTRA{iter}.allMaps(1:(length(runData_EXTRA{iter}.allMaps)/(numStepsToExtract-1)-1):end);
        % end
    end
end

disp('Data loading complete');

%%

fdataRunData = dataRunData;
%fdataRunData(cellfun(@isempty,{dataRunData(:,:,:,:).runData})) = [];
%fdataRunData = reshape(fdataRunData,length(methods),length(mutationRates),length(angles),[]);

disp('Data removal complete');
%%
clear exitCorrectness exitInCorrectness mCor mInCor mCor25 mCor75 mInCor25 mInCor75 userSelDrift
countMaps = zeros(2,1);
for methodID=1:length(methods)
    disp(['Running for method ' int2str(methodID) '/' int2str(length(methods))]);
    for mutRateID=1:length(mutationRates)
    disp(['Running for mutation rate ' int2str(mutRateID) '/' int2str(length(mutationRates))]);
        
        for angle=1:numA
            for exitID=[1:3]
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
    'mutationRates', 'methods', 'angles');

disp('Analysis saved');

save([folder '/' 'selectedMaps.mat'],'dataRunData','-v7.3');
% save([folder '/' 'selectedMaps_allMaps.mat'],'data','-v7.3');

disp('Selected maps saved');

end

