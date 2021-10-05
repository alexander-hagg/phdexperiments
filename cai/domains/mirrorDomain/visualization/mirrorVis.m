function hm = mirrorVis(FV,ffdP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
colors = jet(27);

controlPtsX = 0:1/(ffdP.nDimX-1):1;
controlPtsY = 0:1/(ffdP.nDimY-1):1;
controlPtsZ = 0:1/(ffdP.nDimZ-1):1;
controlPts = combvec(controlPtsX,controlPtsY,controlPtsZ);

allDeforms = permute(ffdP.allDeforms,[2 3 4 5 1]);
deforms = {};
for i=1:3
    for j=1:3
        for h=1:3
            deforms{end+1} = (allDeforms(:,i,j,h,:));
        end
    end
end

controlPtsDelta = cell2mat(deforms);
movedControlpoints = controlPts + controlPtsDelta;

controlPts = controlPts.*repmat((ffdP.maxMeshPoint-ffdP.minMeshPoint)',1,size(controlPts,2)) + repmat(ffdP.minMeshPoint',1,size(controlPts,2));
controlPts = controlPts' * ffdP.rotMat';
controlPts = mapminmax('reverse', controlPts',ffdP.normalizationFactors);

movedControlpoints = movedControlpoints.*repmat((ffdP.maxMeshPoint-ffdP.minMeshPoint)',1,size(movedControlpoints,2)) + repmat(ffdP.minMeshPoint',1,size(movedControlpoints,2));
movedControlpoints = movedControlpoints' * ffdP.rotMat';
movedControlpoints = mapminmax('reverse', movedControlpoints',ffdP.normalizationFactors);

hm = plotmesh(FV.vertices',FV.faces,'facecolor','g');hold on;
% Plot control points
scatter3(controlPts(1,:),controlPts(2,:),controlPts(3,:),32,colors,'filled');hold on;
scatter3(movedControlpoints(1,:),movedControlpoints(2,:),movedControlpoints(3,:),64,colors);
end
