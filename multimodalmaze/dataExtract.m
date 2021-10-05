%function dataExtract(folder,numHidden,numIterations,thetaValues,selCriterion,selValues,penalties)
function dataExtract(folders, numStepsToExtract, recursion)
%%DATAEXTRACT extract partial experimental data, recursion is optional
% 
% folders               - contains a char array of folder names
% numStepsToExtract     - Number of time steps per experiment blob 
%                           if == 1, extract last step
%                           if == 2, extract first and last step
%                           if == n, extract n steps, including first and
%                           last
% recursion             - recurse over folders

% Collect folders recursively
if recursion
    disp('TODO: Implement recursion')
end

for fr=1:length(folders)

end
            name = [int2str(numHidden) '_THETA_' num2str(thetaValues(th)) '_SELCRITID_' ...
                int2str(selCriterion) '_SELVAL_' num2str(selValues(sv)) '_iteration_' int2str(numIterations)];
            disp(name)
            loadedData = load([folder '/out_' name '.mat']); runData = loadedData.runData;
            loadedData = load([folder '/out_' name '_allMaps.mat']); runData_EXTRA = loadedData.runData_EXTRA;
            
            firstGoals{th,sv} = [];
            for runID=2:length(runData)
                for mapID=1:length(runData_EXTRA{runID}.allMaps)
                    map = runData_EXTRA{runID}.allMaps{mapID};
                    firstGoal = map.firstGoal; firstGoal = firstGoal(~isnan(firstGoal));
                    firstGoals{th,sv} = [firstGoals{th,sv};[sum(firstGoal==1),sum(firstGoal==2),sum(firstGoal==3)]];
                end
            end
            save([folder '/data_VALUES_SELCRITID_' int2str(selCriterion) '.mat'],'firstGoals','-v7.3');
            
            % Extract allMaps partially to allow importing on client
            for runID=1:length(runData)
                d{runID} = runData{runID}.d;
                p{runID} = runData{runID}.p;
                if runID==1
                    mapIDs{runID} = [1:2000:length(runData_EXTRA{runID}.allMaps) length(runData_EXTRA{runID}.allMaps)];
                else
                    mapIDs{runID} = [1:200:length(runData_EXTRA{runID}.allMaps) length(runData_EXTRA{runID}.allMaps)];
                end
                for mapID=1:length(mapIDs{runID})
                    allMaps{sv,runID,mapID} = runData_EXTRA{runID}.allMaps{mapIDs{runID}(mapID)};
                end
            end
            save([folder '/data_ALLMAPS_SELCRITID_' int2str(selCriterion) '_THETA_' num2str(thetaValues(th)) '.mat'],'allMaps','mapIDs','d','p','-v7.3');
            
end

end

