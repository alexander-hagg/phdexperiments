function [pltHandle,classification] = viewClasses(samples, classification, varargin)
%viewClasses -
%
% Syntax:  viewClasses(samples,classification,varargin)
%
% Inputs:
%   classStruct - struct  -             
%                           prototypeID: [10 19 9 23 13 4 4 2 10 18 11 5 5 7 3 3 11 3 5 1 2 2 2 1 2]
%                           classMembers: {1?25 cell}
%                           classMembersReduced: {1?25 cell}
%                           prototypes: [25?16 double]
%                           
%
% Outputs:
%   imageHandle - handle of resulting map image
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018

%------------- BEGIN CODE --------------
if nargin > 2; figHandle = varargin{1}; else; figure;figHandle=axes;end

if isempty(classification) 
    classification = extractClasses(samples);
end
    
inds        = classification.simX;
numInds     = size(classification.simX,1);

colors = flipud(colorcube(max(classification.labels))); 
selColors = classification.labels;

%%
pltHandle = scatter(figHandle,inds(:,1), inds(:,2), 8, colors(selColors,:), 'filled');

[~,ids] = unique(classification.labels);

end

