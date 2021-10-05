function [bestIndex, bestBin] = getBestPerCell(fitness, d, edges, features)
%getBestPerCell - Returns index of best individual in each cell
%
% Syntax:  [bestIndex, bestBin] = getBestPerCell(samples,fitness,p)
%
% Inputs:
%   samples - [NXM] - sample genome
%   fitness - [NX1] - fitness values
%   d       - Domain definition struct
%   edges   - {[1 X FeatureRes(1)+1], [1 X FeatureRes(2)+1]}
%           - bin edges of feature space
%
% Outputs:
%   bestIndex - [1XN] - index of best individuals per cell
%
%
% Example:
%
% Other m-files required: d.categorize
%
% See also: createMap

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 27-Jan-2016

%------------- BEGIN CODE --------------

% Get Features of all samples
features(features<0) = 0;
features(features>1) = 1;
for iDim = 1:size(features,2)
    bin(:,iDim) = discretize(features(:,iDim),edges{iDim});
end
[sortedByFeatureAndFitness, indxSortOne] = sortrows([bin fitness],'descend');
[~, indxSortTwo] = unique(sortedByFeatureAndFitness(:,[1 2]),'rows');

bestIndex = indxSortOne(indxSortTwo);
bestBin = bin(bestIndex,:);

%------------- END OF CODE --------------