function folderList = getFoldersRecursive(folder)
%GETFOLDERSRECURSIVE Get folder list recursively
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2019; Last revision: 09-Jan-2019
%
%------------- BEGIN CODE --------------

folderList = [];
fsPaths = dir(folder);
fsPaths(1:2) = [];
for frI=1:length(fsPaths)
    if fsPaths(frI).isdir
        fname = [fsPaths(frI).folder '/' fsPaths(frI).name];
        folderList = [folderList, getFoldersRecursive(fname)];
    else
        folderList = convertCharsToStrings(fsPaths(frI).folder);
    end
end
end

