function parVarMap = getDiversityMap(genesFromMap, varargin)
%GETDIVERSITYMAP Summary of this function goes here
%   Detailed explanation goes here
windowLength = 3; if nargin > 1; windowLength = varargin{1}; end
featureRes(1) = size(genesFromMap,1);featureRes(2) = size(genesFromMap,2);

parVarMap = nan(size(genesFromMap,1),size(genesFromMap,2));
for xx=1:featureRes(1)
    bX = [max(1,xx-windowLength), min(featureRes(1),xx+windowLength)];
    for yy=1:featureRes(2)
        bY = [max(1,yy-windowLength), min(featureRes(1),yy+windowLength)];
        theseGenes = genesFromMap(bX(1):bX(2),bY(1):bY(2),:);
        distances = pdist2(reshape(theseGenes,size(theseGenes,1)*size(theseGenes,2),[]),reshape(theseGenes,size(theseGenes,1)*size(theseGenes,2),[]));
        distances = triu(distances);
        distances = distances(:);
        distances = distances(distances~=0);
        parVarMap(xx,yy) = nanstd(distances);
    end
end
end

