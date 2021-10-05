function dataExtractMultiring(folder, numStepsToExtract, recursion)
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

addpath(genpath('.'));
eval(['delete ' folder '/' 'selectedMaps.mat']);

% Collect folders recursively
if recursion; folders = getFoldersRecursive(folder); end

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
        fname = [folders{fr} '/' fsPaths(fi).name];
        
        angle = regexp(fsPaths(fi).name,'\d*','Match'); angle = str2num(angle{1});
        exitIDs = strfind(folders{fr},'exit');
        if isempty(exitIDs)
            exitNr(1) = 1; exitNr(2) = 1;
        else
            exitNr(1) = str2num(folders{fr}(exitIDs(1)+4)) + 1;
            if length(exitIDs) == 1
                exitNr(2) = 1;
            else
                exitNr(2) = str2num(folders{fr}(exitIDs(2)+4)) + 1;
            end
        end
        
        data(exitNr(1),exitNr(2),angle).fname = fname;
        data(exitNr(1),exitNr(2),angle).description.ringID = exitNr(1)-1;
        data(exitNr(1),exitNr(2),angle).description.exitID = exitNr(2)-1;
        data(exitNr(1),exitNr(2),angle).description.angle = angle;
        load(fname,'runData_EXTRA');
        for iter=1:length(runData_EXTRA)
            if isempty(runData_EXTRA{iter}); continue; end
            if numStepsToExtract == 1
                data(exitNr(1),exitNr(2),angle).maps(iter,:) = runData_EXTRA{iter}.allMaps(end);
            else
                data(exitNr(1),exitNr(2),angle).maps(iter,:) = runData_EXTRA{iter}.allMaps(1:length(runData_EXTRA{iter}.allMaps)/(numStepsToExtract-1)-1:end);
            end
        end
    end
end

% data = data(~cellfun(@isempty,{data.fname}));

save([folder '/' 'selectedMaps.mat'],'data','-v7.3');
disp('SAVED and DONE');

end

