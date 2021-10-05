function vTrajectories(trajectories, values, maze, drawMaze, varargin)
%visTrajectories
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018
%
%------------- BEGIN CODE --------------
showVal = 3; if nargin>4;showVal=varargin{1};end
hold off;
if drawMaze
    mapImgInput = importdata(['mazeDomain/simulator/worlds/' maze '.pbm']);
    mapImgInput = ~mapImgInput;
    image(mapImgInput);
    colormap([1,1,1;0,0,0]);
    hold on;
end

if showVal == 1
    colorScaling = values>1;
    trajColor = ((1-colorScaling).*[0 1 0]')' + (colorScaling.*[1 0 0]')';
elseif showVal == 2 || showVal == 3
    colorScaling = values;
    clrs = [parula(4)];
    trajColor = clrs(colorScaling+1,:);
else
    trajColor = [0 0 0];
end

for ii=1:size(trajectories,1)
    invalidTrajSegment = trajectories(ii,:,1)~=0;
    
    if size(trajColor,1) < size(trajectories,1)
        plot(trajectories(ii,invalidTrajSegment,1)',trajectories(ii,invalidTrajSegment,2)','Color',trajColor,'LineWidth',1);
        scatter(trajectories(ii,end,1)',trajectories(ii,end,2)',64,trajColor,'filled');
    else
        plot(trajectories(ii,invalidTrajSegment,1)',trajectories(ii,invalidTrajSegment,2)','Color',trajColor(ii,:),'LineWidth',1);
        scatter(trajectories(ii,end,1)',trajectories(ii,end,2)',64,trajColor(ii,:),'filled');
    end
end


axis([0 400 0 400]);
grid on;
ax = gca;ax.XTickLabel = [];ax.YTickLabel = [];

if ~drawMaze; hold on; end;

end

