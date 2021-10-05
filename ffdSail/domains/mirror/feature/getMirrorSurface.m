function area = getMirrorSurface(mirror, d)
%getMirrorSurface - surface area of reflective area of mirror
%
% Syntax:  [area] = getMirrorSurface(mirror, d)
%
% Inputs:
%    mirror - [Npoints X 3] - X,Y,Z coordinates of each changeable vertex in design
%    d      - domain
% Outputs:
%    area - [scalar] - Area of mirror 
%
% Example:
%    d = mirror_Domain; % for point IDs
%    FV = mirror_ffd_Express(0.5 + zeros(1,41), d.FfdP);
%    area = getMirrorSurface(FV.vertices,d)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: alexander.hagg@h-brs.de
% Dec 2017; Last revision: 12-Dec-2017

%% ------------- BEGIN CODE --------------

mirrorSubmesh = mirror(:,d.features.mirror.ids);
mirrorSubmesh = mirrorSubmesh - mean(mirrorSubmesh')';
mirrorSubmesh = (mirrorSubmesh'*d.FfdP.rotMat)';
%%
area = polyarea(mirrorSubmesh(2,:)',mirrorSubmesh(3,:)');
%------------- END OF CODE --------------











