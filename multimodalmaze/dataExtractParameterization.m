function dataExtractParameterization(folder, numStepsToExtract, recursion)
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
        fname = [folders{fr} '/' fsPaths(fi).name]
        
        angle = regexp(fsPaths(fi).name,'\d*','Match'); angle = str2num(angle{1});
        pweightIDs = strfind(folders{fr},'pweight');
        if ~isempty(pweightIDs)
            pweightNR = str2num(folders{fr}(pweightIDs(1)+7:end));
        else
            continue;
        end

        data(pweightNR*1000+1,angle).fname = fname;
        data(pweightNR*1000+1,angle).description.penaltyweight = pweightNR;
        data(pweightNR*1000+1,angle).description.angle = angle;
        load([fsPaths(fi).folder '/' int2str(angle) '.mat'],'runData');
        data(pweightNR*1000+1,angle).runData = runData;
        
        load(fname,'runData_EXTRA');
        for iter=1:length(runData_EXTRA)
            if isempty(runData_EXTRA{iter}); continue; end
            data_EXTRA(pweightNR*1000+1,angle).maps(iter,:) = runData_EXTRA{iter}.allMaps(1:(length(runData_EXTRA{iter}.allMaps)/(numStepsToExtract-1)-1):end);
        end        
        
    end
end

save([folder '/' 'selectedMaps.mat'],'data','-v7.3');
save([folder '/' 'selectedMaps_allMaps.mat'],'data_EXTRA','-v7.3');

disp('SAVED and DONE');

end

