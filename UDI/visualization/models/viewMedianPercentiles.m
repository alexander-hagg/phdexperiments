function hndl = viewMedianPercentiles( inMED, inPRC, varargin )
%VIEWMEDIANPERCENTILES Summary of this function goes here
%   Detailed explanation goes here

% Parse input arguments
parse = inputParser;
parse.addRequired('inMED');
parse.addRequired('inPRC');
parse.addOptional('tickLen',4,@isnumeric);
parse.addOptional('tickRelDis',0.02,@isnumeric);
parse.addOptional('groupDist',8,@isnumeric);
parse.addOptional('tickLW',2,@isnumeric);
parse.addOptional('tickVertLW',2,@isnumeric);
parse.addOptional('medianLW',2,@isnumeric);
parse.addOptional('groupNames',[]);
parse.addOptional('groupSelected',[]);
parse.addOptional('subgroups',[]);
parse.addOptional('prcNames',[]);
parse.addOptional('colorScheme', 'hsv', @ischar);

parse.parse(inMED, inPRC, varargin{:});
tickLength          = parse.Results.tickLen;
tickRelDis          = parse.Results.tickRelDis;
groupDist           = parse.Results.groupDist;
tickLineWidth       = parse.Results.tickLW;
tickVertLineWidth   = parse.Results.tickVertLW;
medianLineWidth     = parse.Results.medianLW;
groupNames          = parse.Results.groupNames;
groupSelected       = parse.Results.groupSelected;
subgroups           = parse.Results.subgroups;
prcNames            = parse.Results.prcNames;
colorScheme         = parse.Results.colorScheme;

% Select data
if ~isempty(groupSelected)
    inMED               = inMED(:,groupSelected);
    inPRC               = inPRC(:,:,groupSelected);
    groupNames = groupNames(groupSelected);
end

% Use some sane color
if isempty(subgroups)
    clrs = feval(colorScheme, size(inMED,1)+2);
    clrs = clrs(2:size(inMED,1)+1,:);
else
    clrs = feval(colorScheme, size(subgroups,2)+1);    
    clrs = clrs(1:size(inMED,1),:);
end

%% Calculate coordinates for vertical ticks
ycoords = inPRC;
xc = groupNames;
if isempty(subgroups)
    disp('Error: define subgroups');
else
    delta = 0:tickRelDis:tickRelDis*(length(subgroups)-1);
end

xx(1,:,:) = xc + (xc'*delta)';
xcoords = repmat(xx,size(inPRC,1),1,1);

%% Fix color order
co = get(gca,'ColorOrder');
set(gca, 'ColorOrder', repelem(clrs,size(inPRC,3),1), 'NextPlot', 'replacechildren');

% Plot vertical ticks
for i=1:size(xcoords,2)
    for j=1:size(xcoords,3)
        plot(squeeze(xcoords(:,i,j)),squeeze(ycoords(:,i,j)),'-','LineWidth',tickVertLineWidth);
        hold on;
    end    
end

%%
% Calculate coordinates for horizontal ticks and plot
% ticksYU = repelem(ycoords(1,:),2,1);
% ticksXU = [xcoords(1,:)-tickLength*ones(size(xc)) ;xcoords(2,:)+tickLength*ones(size(xc)) ];
% plot(ticksXU,ticksYU,'LineWidth',tickLineWidth);
% 
% ticksYD = repelem(ycoords(2,:),2,1);
% plot(ticksXU,ticksYD,'LineWidth',tickLineWidth);

%%
for i=1:size(xcoords,2)
    hndl(i,1) = plot(squeeze(xcoords(1,i,:)),inMED(i,:), '-', 'Color', clrs(i,:), 'LineWidth',medianLineWidth);
    scatter(squeeze(xcoords(1,i,:)),inMED(i,:), 64, clrs(i,:), 'filled');
end

grid on;

%max([ticksYU(:);ticksYD(:)]);
%min([ticksYU(:);ticksYD(:)]);
ax = gca;
ax.XTick = groupNames;
if ~isempty(groupNames);ax.XTickLabel = string(groupNames);end

end

