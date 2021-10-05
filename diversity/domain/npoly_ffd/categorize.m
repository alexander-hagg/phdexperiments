function features = categorize(phenotypes, d)
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

for i=1:length(phenotypes)
    pgon = phenotypes{i};
    
    % Feature 1: Area
    features(i,1) = area(pgon);
    
    % Feature 2: Perimeter
    features(i,2) = perimeter(pgon);
    
    %vertices = pgon.Vertices(all(~isnan(pgon.Vertices)'),:);    
    %vertices = unique(vertices,'rows','stable');
    %y = interppolygon([vertices],100,'linear');
    %vertexDistances = mypdist2(y,y);
    %vertexDistances(logical(eye(size(vertexDistances,1)))) = nan;
    %features(i,3) = nanmax(vertexDistances(:));
    %features(i,4) = nanmin(vertexDistances(:));
end

features(:,1) = (features(:,1)-d.featureMin(1))./(d.featureMax(1)-d.featureMin(1));
features(:,2) = (features(:,2)-d.featureMin(2))./(d.featureMax(2)-d.featureMin(2));
%features(:,3) = (features(:,3)-d.featureMin(3))./(d.featureMax(3)-d.featureMin(3));
%features(:,4) = (features(:,4)-d.featureMin(4))./(d.featureMax(4)-d.featureMin(4));

% Add feature "random" for debugging purposes
%features(:,5) = rand(size(features(:,1)));

features(features>1) = 1; features(features<0) = 0;

end

%------------- END CODE -------------- 