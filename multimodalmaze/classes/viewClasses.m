function [pltHandle] = viewClasses(classification,showLabels,varargin)
%viewClasses -
%
% Syntax:  viewClasses(classStruct,d,varargin)
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

inds        = classification.simX;
numInds     = size(classification.simX,1);

if nargin > 2; selection=varargin{1};end
if exist('selection','var')
    disp('Showing selected classes');
    colors = [1 0 0;0 1 0];
    selColors = ismember(classification.labels,selection)+1;
else
    colors = flipud(colorcube(max(classification.labels))); 
    selColors = classification.labels;
end

%%
pltHandle = scatter(inds(:,1), inds(:,2), 32, colors(selColors,:), 'filled');

[~,ids] = unique(classification.labels);

if showLabels
    text(inds(ids,1)+0.1,inds(ids,2)+0.1,[int2str(classification.labels(ids))],'FontSize',16);
end

end

