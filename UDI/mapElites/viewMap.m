function [figHandle, imageHandle, cHandle] = viewMap(mapMatrix, d, varargin)
%computeFitness - Computes fitness with penalties from drag, lift, area
%
% Syntax:  viewMap(predMap.fitness, d)
%
% Inputs:
%   mapMatrix   - [RXC]  - scalar value in each bin (e.g. fitness)
%   d           - struct - Domain definition
%
% Outputs:
%   figHandle   - handle of resulting figure
%   imageHandle - handle of resulting map image
%
%
% Example:
%    p = sail;
%    d = af_Domain;
%    output = sail(d,p);
%    d.featureRes = [50 50];
%    predMap = createPredictionMap(output.model,p,d);
%    viewMap(predMap.fitness,d)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: createMap, updateMap, createPredictionMap

% Author: Adam Gaier, Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 19-Jun-2019

%------------- BEGIN CODE --------------
mapRes = size(mapMatrix);
for i=1:length(mapRes)
    edges{i} = linspace(0,1,mapRes(i)+1); %#ok<AGROW>
end

yOffset = [0.5 -0.0 0];
imgHandle = imagesc(flipud(rot90(mapMatrix))); fitPlot = gca;
if nargin > 3
    if strcmp(varargin{1},'flip')
        imgHandle = imagesc(fliplr(rot90(rot90(mapMatrix)))); fitPlot = gca;
    end
end

set(imgHandle,'AlphaData',~isnan(imgHandle.CData)*1)
xlab = xlabel([d.featureLabels{1} '\rightarrow']);
ylab = ylabel(['\leftarrow' d.featureLabels{2} ]);
%set(ylab,'Rotation',0,'Position',get(ylab,'Position')-yOffset)


xticklabels = {}; yticklabels = {};
set(fitPlot,...
    'XTickLabel',xticklabels,...
    'XTick', linspace(0.5,mapRes(2)+0.5,mapRes(2)+1), ...
    'YTickLabel',yticklabels,...
    'YTick', linspace(0.5,mapRes(1)+0.5,mapRes(1)+1),...
    'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-',...
    'xcolor', 'k', 'ycolor', 'k'...
    )

cHandle = colorbar;
axis square
figHandle = fitPlot; imageHandle = imgHandle;

%------------- END OF CODE --------------