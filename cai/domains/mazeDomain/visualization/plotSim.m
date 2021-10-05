function plotSim(trajectory,d,varargin)
%PLOTSIM Summary of this function goes here
%%   Detailed explanation goes here
plotmap = true; if nargin > 2; plotmap = varargin{1};end
clr = [0.7 0.7 0.7]; if nargin > 3; clr = varargin{2};end
showPath = true; if nargin > 4; showPath = varargin{3};end
if plotmap
    img = imread('mazeDomain/simulator/worlds/mediumRound.pbm');
    imagesc(img); colormap([0,0,0;1,1,1]);
    hold on;
    scatter(d.goalLocations(:,1),d.goalLocations(:,2),256,[1 0 0;0 1 0; 0 0 1],'x','LineWidth',4);
    axis([-100 500 -100 500]);
end

if ndims(trajectory)==2
    if showPath; pHandle = plot(trajectory(:,1),trajectory(:,2),'-','LineWidth',1,'Color',[0 0 0]);end
    scatter(trajectory(end,1),trajectory(end,2),32);
else
    if showPath; pHandle = plot(trajectory(:,:,1)',trajectory(:,:,2)','-','LineWidth',1,'Color',[0 0 0]);end
    scatter(trajectory(:,end,1)',trajectory(:,end,2)',32,clr,'filled');
end
hold on;

%%
end
