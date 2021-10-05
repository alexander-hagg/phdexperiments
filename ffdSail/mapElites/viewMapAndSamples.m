function [figHandle, imageHandle, cHandle] = viewMapAndSamples(mapMatrix, d, varargin)
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

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 20-Aug-2017

%------------- BEGIN CODE --------------

parse = inputParser;
parse.addRequired('mapMatrix');
parse.addRequired('d');
parse.addOptional('samples', []);
parse.addOptional('sampleSymbolSize', 10);

parse.parse(mapMatrix,d,varargin{:});
mapMatrix           = parse.Results.mapMatrix;
d                   = parse.Results.d;
samples             = parse.Results.samples;
sampleSymbolSize    = parse.Results.sampleSymbolSize;

mapRes = size(mapMatrix);
for i=1:length(mapRes)
    edges{i} = linspace(0,1,mapRes(i)+1); %#ok<AGROW>
end
%%
yOffset = [0.5 -0.0 0];
imgHandle = imagesc(flipud(rot90(mapMatrix))); fitPlot = gca;
hold on;

if ~isempty(samples)
    %%
    for i=1:size(samples,1)
        [II,JJ] = ind2sub(size(mapMatrix),samples(i));
        coord(i,:) = [size(mapMatrix,1) + 1 - II,JJ];
    end
    indices = sub2ind(size(mapMatrix),coord(:,1),coord(:,2));
    totals = accumarray(indices(:), 1);
    if length(totals) < size(mapMatrix,1)*size(mapMatrix,2)
        totals(length(totals)+1:size(mapMatrix,1)*size(mapMatrix,2)) = 0;
    end
    
    [XX, YY] = ind2sub(size(mapMatrix),1:625);
    XX(totals==0) = [];YY(totals==0) = [];
    
    scatter(XX(:), YY(:), sampleSymbolSize*totals(totals~=0), 'r','filled');
    
end
%%
set(imgHandle,'AlphaData',~isnan(imgHandle.CData)*1)
xlab = xlabel([d.featureLabels{1} '\rightarrow']);
ylab = ylabel(['\downarrow' d.featureLabels{2} ]);
set(ylab,'Rotation',0,'Position',get(ylab,'Position')-yOffset)



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
set(fitPlot,...
    'XTickLabel',xticklabels,...
    'XTick', linspace(0.5,d.featureRes(2)+0.5,d.featureRes(2)+1), ...
    'YTickLabel',yticklabels,...
    'YTick', linspace(0.5,d.featureRes(1)+0.5,d.featureRes(1)+1),...
    'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-',...
    'xcolor', 'k', 'ycolor', 'k'...
    )

cHandle = colorbar;
axis square
figHandle = fitPlot; imageHandle = imgHandle;

%------------- END OF CODE --------------