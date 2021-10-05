function [features,nonNormalizedFeatures] = categorize(phenotypes, d)
%categorize - extracts all features of the npoly domain
%
% Syntax:  features = categorize(samples, phenotypes, d)
%
% Inputs:
%    samples        - [NxD] - N samples with dimensionality D. Sometimes
%                             used to circumvent having to instantiate full
%                             phenotype
%    phenotypes     - [cell array] - phenotypes (expressed samples)
%    d              - struct  - Domain description struct
% 
%
% Outputs:
%    features       - [NxF] - F features for N samples
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE -------------- 

if isempty(phenotypes)
    features = [];
    return;
end

for i=1:length(phenotypes)
    pgon = phenotypes{i};
    pgon = simplify(pgon);
    
    % Feature 1: Area
    nonNormalizedFeatures(i,1) = area(pgon);
    
    % Feature 2: Perimeter
    nonNormalizedFeatures(i,2) = perimeter(pgon);
end

features(:,1) = (nonNormalizedFeatures(:,1)-d.featureMin(1))./(d.featureMax(1)-d.featureMin(1));
features(:,2) = (nonNormalizedFeatures(:,2)-d.featureMin(2))./(d.featureMax(2)-d.featureMin(2));

features(features>1) = 1; features(features<0) = 0;

end

%------------- END CODE -------------- 
