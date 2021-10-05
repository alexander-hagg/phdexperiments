function curvature = getTotalCurvature(mirror, d)
%getTotalCurvature - Returns total curvature of a set of 2D lines
%
% Syntax:  [totalCurvature] = getTotalCurvature(mirror, d)
%
% Inputs:
%    mirror     - [Npoints X 3] - X,Y,Z coordinates of each point in design
%    d - domain      - Ids of points in line for curvature
%
% Outputs:
%    meanCurvature - [scalar] - Mean curvature of all lines in 2D
%
% Example:
%   d = mirror_Domain; % for point IDs
%   FV = mirror_ffd_Express(0.5 + zeros(1,41), d.FfdP)
%   meanCurvature = getTotalCurvature(FV.vertices, d.curvSecIds)
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

% Get length relative curvature and length of each line
p1 = mirror(:,d.features.curvature.ids(1));
p2 = mirror(:,d.features.curvature.ids(12));
theta           = -atan((p1(1)-p2(1))/(p1(2)-p2(2)));
rotMat          = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
mirror = (mirror'*rotMat)';
line = mirror(:,d.features.curvature.ids)';
lineCurv = LineCurvature2D(line(:,[2 3]));
curvature = sum(abs(lineCurv(:)));

%% Visualize mirror, curvature line and extremes (which determined angle)
% figure(7);hold off;
% scatter3(mirror(1,:),mirror(2,:),mirror(3,:)); hold on;
% view(0,90)
% plot3(line(:,1),line(:,2),line(:,3),'k-','LineWidth',4);
% id = 1
% scatter3(line(id,1),line(id,2),line(id,3),256,'filled');
% 
% id = 12
% scatter3(line(id,1),line(id,2),line(id,3),256,'filled');

end
%------------- END OF CODE --------------











