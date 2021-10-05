function mapLinIndx = getBins( inputs, coordFunc, d, edges )
%GETBINS Summary of this function goes here
%   Detailed explanation goes here
sampleCoords = coordFunc(inputs);
for iDim = 1:d.nDims
    sampleBin(:,iDim) = discretize(sampleCoords(:,iDim), edges{iDim});
end
mapLinIndx = sub2ind(d.featureRes,sampleBin(:,1),sampleBin(:,2));

end

