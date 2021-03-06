function [axHandle, imageHandle, cHandle] = viewMap(mapMatrix, d, varargin)
%
% Syntax:  [axHandle, imageHandle, cHandle] = viewMap(mapMatrix, d, varargin)
%
% Inputs:
%   mapMatrix   - [RXC]  - scalar value in each bin (e.g. fitness)
%   d           - struct - Domain definition
%
% Outputs:
%   axHandle   - handle of resulting figure
%   imageHandle - handle of resulting map image
%
%
% Example:
%    p = sail;
%    d = af_Domain;
%    output = sail(d,p);
%    d.featureResolution = [50 50];
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
% email: adam.gaier@h-brs.de, alexander.hagg@h-brs.de
% Jun 2016; Last revision: 15-Aug-2019

%------------- BEGIN CODE --------------
if nargin > 2; mapValue = varargin{1}; else; mapValue = 'fitness'; end
if nargin > 3; axHandle = varargin{2}; else; axHandle = gca;end
hold(axHandle,'off');
mapMatrix = eval(['mapMatrix.' mapValue]);
mapMatrix = squeeze(mapMatrix);
mapRes = size(mapMatrix);
for i=1:length(mapRes)
    edges{i} = linspace(0,1,mapRes(i)+1); %#ok<AGROW>
end

imgHandle = imagesc(axHandle,fliplr(rot90(rot90(rot90(mapMatrix)))));

set(imgHandle,'AlphaData',~isnan(imgHandle.CData)*1)
xlab = xlabel(axHandle,[d.featureLabels{d.selectedFeatures(1)} '\rightarrow']);
ylab = ylabel(axHandle,['\leftarrow' d.featureLabels{d.selectedFeatures(2)} ]);


xticklabels = num2str(edges{2}',2);
if length(xticklabels)>10
    keep = false(1,length(xticklabels));
    keep(1:round(length(xticklabels)/10):length(xticklabels)) = true;
    xticklabels(~keep,:)= ' ';
end

yticklabels = num2str(edges{1}(end:-1:1)',2);
if length(yticklabels)>10
    keep = false(1,length(yticklabels));
    keep(1:round(length(yticklabels)/10):length(yticklabels)) = true;
    yticklabels(~keep,:)= ' ';
end
xticklabels = {}; yticklabels = {};
set(axHandle,...
    'XTickLabel',xticklabels,...
    'XTick', linspace(0.5,size(mapMatrix,2)+0.5,size(mapMatrix,2)+1), ...
    'YTickLabel',yticklabels,...
    'YTick', linspace(0.5,size(mapMatrix,1)+0.5,size(mapMatrix,1)+1),...
    'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-',...
    'xcolor', 'k', 'ycolor', 'k'...
    )
cmap = parula(33); cmap(end,:) = [];
colormap(axHandle,cmap);
cHandle = colorbar(axHandle);
axis(axHandle,'square');
axis(axHandle,[0.5 size(mapMatrix,1)+0.5 0.5 size(mapMatrix,2)+0.5]);

imageHandle = imgHandle;
hold(axHandle,'on')
%features = reshape(mapMatrix.features,mapRes(1)*mapRes(2),[]);
%scatter(axHandle,0.5+features(:,d.selectedFeatures(1))*mapRes(1),0.5+features(:,d.selectedFeatures(2))*mapRes(2),8,[0 0 0],'filled')
%------------- END OF CODE --------------