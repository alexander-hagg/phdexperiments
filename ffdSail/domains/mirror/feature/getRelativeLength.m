function [length, maxVertex, minVertex, maxVertexID, minVertexID] = getLength(mirror, d)
%getRelativeLength - Returns relative length of mirror (x/y)
%
% Syntax:  [relativeLength] = getLength(mirror, d)
%
% Inputs:
%    mirrorSubmesh - [Npoints X 3] - X,Y,Z coordinates of each changeable vertex in design
%    d - domain
% Outputs:
%    length - [scalar] - Relative length of mirrorSubmesh (x/y)
%
% Example:
%    d = mirror_Domain; % for point IDs
%    FV = mirror_ffd_Express(0.5 + zeros(1,41), d.FfdP);
%    length = getLength(FV.vertices,d)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: alexander.hagg@h-brs.de
% Dec 2017; Last revision: 12-Dec-2017

%%------------- BEGIN CODE --------------
mirror = mirror - mean(mirror')';
mirror = mirror'*d.FfdP.rotMat;
mirrorSubmesh = mirror(d.features.subMeshIds,:);
[maxVertex,maxVertexID] = max(mirrorSubmesh(:,1));
[minVertex,minVertexID] = min(mirrorSubmesh(:,1));
[~,maxVertexID] = intersect(mirror(:,2)',mirrorSubmesh(maxVertexID,2));
[~,minVertexID] = intersect(mirror(:,2)',mirrorSubmesh(minVertexID,2));
length = (maxVertex - minVertex);
%------------- END OF CODE --------------











