function showHistory(constraints,classification,d,varargin)
%showHistory - show selection history
%
% Syntax:  showHistory(constraints,classification,d,handle)
%
% Inputs:
%    classification
%    selection
%    type
%    method
%
% Outputs:
%    constraints
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Aug 2019; Last revision: 16-Aug-2019
%
%------------- BEGIN CODE --------------
if nargin > 3; handle = varargin{1}; else; figure; handle = axes; end
numCols = length(constraints);
for i=1:numCols
    numPrototypes = length(constraints{i}.selectedClasses);
    placement = [];
    for j=1:numPrototypes
        placement(j,:) = [(i-1) (j-1)]*d.spacer;
    end
    showPhenotype(handle,d,classification{i}.protoX(constraints{i}.selectedClasses,:),placement);
end
axis(handle,'equal');
title(handle,'Selection History');
end

