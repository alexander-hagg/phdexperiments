function [figHandle, imageHandle, cHandle] = viewMap(map, d, varargin)
%computeFitness - Computes fitness with penalties from drag, lift, area
%
% Syntax:  [figHandle, imageHandle, cHandle] = viewMap(mapMatrix, d, varargin)
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: createMap, updateMap, createPredictionMap

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jul 2019; Last revision: 15-Aug-2019

%------------- BEGIN CODE --------------
if nargin > 2; figHandle = varargin{1}; else; figHandle = figure;end

elites = map.features(:,d.selectedFeatures);
[elites,ids] = unique(elites,'rows');
fitness = map.fitness(ids);
hold(figHandle,'off');
h = imagesc(figHandle,1);delete(h);
[v,c]=voronoin(elites);
imageHandle = voronoi(figHandle,elites(:,1),elites(:,2));
v1 = shiftdim(reshape([imageHandle(2).XData;imageHandle(2).YData],2,3,[]),2); % Arranged one edge per row, one vertex per slice in the third dimension
nUnbounded = sum(cellfun(@(ic)ismember(1,ic),c));
v1Unbounded = v1(end-(nUnbounded-1):end,:,:);
[~,iBounded] = min(pdist2(v,v1Unbounded(:,:,1))); % Index of the bounded vertex
vUnbounded = v1Unbounded(:,:,2); % Displayed coordinate of the unbounded end of the cell edge

l = 0;
maxPatch = 100;
patchesX = nan(maxPatch,size(elites,1));patchesY = nan(maxPatch,size(elites,1));
for s=1:size(elites,1)
    l=l+1;
    if l > length(c); continue;end
    cPatch = c{l}; % List of vertex indices
    vPatch = v(cPatch,:); % Vertex coordinates which may contain Inf
    idx = find(cPatch==1); % Check if cell has unbounded edges
    if idx
        cPatch = circshift(cPatch,-idx); % Move the 1 to the end of the list of vertex indices
        vPatch = [vPatch(1:idx-1,:)
            vUnbounded(iBounded == cPatch(end-1),:)
            vUnbounded(iBounded == cPatch(1),:)
            vPatch(idx+1:end,:)]; % Replace Inf values at idx with coordinates from the unbounded edges that meet the two adjacent finite vertices
    end
    vPatch = padarray(vPatch,[maxPatch-size(vPatch,1) 1],'replicate','post');
    patchesX(1:length(vPatch(:,1)),s) = vPatch(:,1);
    patchesY(1:length(vPatch(:,2)),s) = vPatch(:,2);
end

hold(figHandle,'on');
patch(figHandle,patchesX,patchesY,fitness);
cmap = hot(33); cmap(end,:) = [];
colormap(figHandle,cmap);
cHandle = colorbar(figHandle);
cHandle.Label.String = 'Fitness';
axis(figHandle,[0 1 0 1]);

scatter(figHandle,elites(:,1),elites(:,2),8,[0 0 0],'filled');

xlab = xlabel(figHandle,[d.featureLabels{d.selectedFeatures(1)} '\rightarrow']);
ylab = ylabel(figHandle,['\leftarrow' d.featureLabels{d.selectedFeatures(2)} ]);

%------------- END OF CODE --------------