function hm = mirrorVisPaper(FV,ffdP, d, rotation, viewRotation, showCtlPts, showMovedCtlPts, annotate, FVRED, showFeature1, showFeature2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4;rotation = 0;end
if nargin < 5;viewRotation = false;end
if nargin < 6;showCtlPts = false;end
if nargin < 7;showMovedCtlPts = false;end
if nargin < 8;annotate = false;end
if nargin < 10;showFeature1 = false;end
if nargin < 11;showFeature2 = false;end
colors = jet(27);

controlPtsX = 0:1/(ffdP.nDimX-1):1;
controlPtsY = 0:1/(ffdP.nDimY-1):1;
controlPtsZ = 0:1/(ffdP.nDimZ-1):1;
controlPts = combvec(controlPtsX,controlPtsY,controlPtsZ);

allDeforms = ffdP.allDeforms;
deforms = {};
for i=1:ffdP.nDimY
    for j=1:ffdP.nDimX
        for h=1:ffdP.nDimZ
            deforms{end+1} = (allDeforms(:,i,j,h,:));
        end
    end
end

%%
if exist('FVRED','var')
    F = FVRED.faces;
    V = FVRED.vertices;    
else
    F = FV.faces;
    V = FV.vertices';    
end

hm = 0;
%k = boundary(V(d.features.curvature.ids,1),V(d.features.curvature.ids,2),V(d.features.curvature.ids,3));
hm = plotmesh(V,F);hold on;
k = d.features.curvature.ids;
plot3(V(k,1),V(k,2),V(k,3),'k:','LineWidth',4);
hold on;

rotMat = [cos(rotation),-sin(rotation),0; ...
          sin(rotation), cos(rotation),0; ...
          0            , 0            ,1];
      
      
VT = V*rotMat;
if viewRotation==0
    k = boundary(VT(:,2),VT(:,3));
elseif viewRotation==90
    k = boundary(VT(:,1),VT(:,2));   
end
plot3(V(k,1),V(k,2),V(k,3),'r-','LineWidth',4);

%hm.EdgeColor = [0.3 0.3 0.3];
%hm.FaceColor = [1 1 1];
%% Plot control points

if showCtlPts && ~showMovedCtlPts
    controlPts = ffdP.controlPts;
    normFac = ffdP.normalizationFactors;
    controlPts = repmat(normFac.xmin,1,size(controlPts,2)) + (controlPts-normFac.ymin).*repmat(normFac.xrange,1,size(controlPts,2))./normFac.yrange;

    scatter3(controlPts(1,:),controlPts(2,:),controlPts(3,:),32,[0 0 0],'filled');
    hold on;
end
% if showMovedCtlPts
%     xxx = [controlPts(1,:);movedControlpoints(1,:)];
%     yyy = [controlPts(2,:);movedControlpoints(2,:)];
%     zzz = [controlPts(3,:);movedControlpoints(3,:)];
%     plot3(xxx,yyy,zzz,'LineWidth',2);hold on;
%     scatter3(xxx(1,:),yyy(1,:),zzz(1,:),'filled');
%     scatter3(xxx(2,:),yyy(2,:),zzz(2,:),'<','filled');
% %    scatter3(movedControlpoints(1,:),movedControlpoints(2,:),movedControlpoints(3,:),48,[1 0 0],'filled');
% end


sel{1} = [1 3 12 10 1];sel{2} = [25 27 36 34 25];sel{3} = [1 25];sel{4} = [3 27];sel{5} = [12 36];sel{6} = [10 34];
if showCtlPts
    for i=1:length(sel)
        plot3(controlPts(1,sel{i}),controlPts(2,sel{i}),controlPts(3,sel{i}),'k','LineWidth',1);hold on;
    end
end
if annotate
    if showCtlPts
        for i=1:size(controlPts,2)
            a = text(controlPts(1,i)-10,controlPts(2,i)-15,controlPts(3,i),int2str(i))
            a.FontSize = 18;
        end
    end
    if showMovedCtlPts
        for i=1:size(controlPts,2)
            a = text(movedControlpoints(1,i)-10,movedControlpoints(2,i)-15,movedControlpoints(3,i),[int2str(i) '"'])
            a.FontSize = 18;
        end
    end
end


if showFeature1
    % Plot curvature lines
    for i=1:length(d.features.curvature.ids)
        pp = FV.vertices(:,d.features.curvature.ids{i});
        h(2) = plot3(pp(1,:),pp(2,:),pp(3,:),'r:','LineWidth',2);
    end
end

if showFeature2
    [~,~, ~, maxVertexID, minVertexID] = getRelativeLength(FV.vertices, d);
    % Plot relative length
    pp = FV.vertices(:,[maxVertexID]);
    pp(3,:) = 1000;
    h(3) = scatter3(pp(1,:),pp(2,:),pp(3,:),192,'w','filled');
    h(3) = scatter3(pp(1,:),pp(2,:),pp(3,:),128,'k','filled');
    h(3) = scatter3(pp(1,:),pp(2,:),pp(3,:),64,'r','filled');
end

end
