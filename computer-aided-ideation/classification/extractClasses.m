function [classification,simspaceCoordinates] = extractClasses(samples,varargin)
%extractClasses Extract classes, determine prototypes based on fitness,
%add extra columns
%
% Syntax:  [classification,simX] = extractClasses(samples,varargin)
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

clusterMethod           = 'dbscan';
if nargin > 3
    clusterMethod = varargin{2};
    numClusters = varargin{3};
end

% Reshape data
if ndims(samples) > 2
    samples = reshape(samples,size(samples,1)*size(samples,2),[]);
    samples = samples(all(~isnan(samples')),:);
end

if nargin > 1 && ~isempty(varargin{1})
    simspaceCoordinates = varargin{1};
else
    simspaceCoordinates             = getSimSpace(samples);
end

if strcmp(clusterMethod,'dbscan')
    coreneighbours      = max(2 * numDims_DR,3); %Rule of thumb
    [~,t_distances]     = knnsearch(simspaceCoordinates,simspaceCoordinates,'K',coreneighbours+1);
    t_distances(:,1)    = [];
    t_distances         = sort(t_distances(:));
    [maxVal ,~]         = getElbow(t_distances);
    epsilon             = maxVal;
    [~,labels,cen] = dbscan(simspaceCoordinates', epsilon, coreneighbours);
    
elseif strcmp(clusterMethod,'kmedoids')
    [labels,cen] = kmedoids(simspaceCoordinates,numClusters);
end

if strcmp(clusterMethod,'dbscan')
    labels = labels + 1; % Get rid of zero label
end

% Prototypes and classes
for iii=1:length(labels)
    ii = labels(iii);
    Xclass = samples(labels==ii,:); simXclass = simspaceCoordinates(labels==ii,:);
    
    id = ismember(simXclass,cen(ii,:),'rows');
    classification.protoX(ii,:) = Xclass(id,:); classification.protoSimX(ii,:) = simXclass(id,:);
end

classification.X = samples;
classification.simX = simspaceCoordinates;
classification.labels = labels;

end

%------------- END CODE --------------
