function simX = getSimSpace(samples)
%getSimSpace - get similarity spaces coordinates with tSNE
%
% Syntax:  simX = getSimSpace(X)
%
% Inputs:
%    samples - coordinates in original space
%
% Outputs:
%    classification
%    simspaceCoordinates
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018
%
%------------- BEGIN CODE --------------

% Reshape data
if ndims(samples) > 2
    samples = reshape(samples,size(samples,1)*size(samples,2),[]);
    samples = samples(all(~isnan(samples')),:);
end

numDims_DR              = 2;
numDims_ORG             = size(samples,2);
numSamples              = size(samples,1);
perplexity              = min(floor(numSamples*0.33),50);
speedQualitytradeOff    = 0.3;

simX = fast_tsne(samples, numDims_DR, numDims_ORG, perplexity, speedQualitytradeOff);

end

%------------- END CODE --------------
