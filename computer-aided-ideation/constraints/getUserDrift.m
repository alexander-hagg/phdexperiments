function [drift,simspacePredictions] = getUserDrift(samples,c)
%getUserDrift - calculate user drift
%
% Syntax:  [drift,simspacePredictions] = getUserDrift(samples,c)
%
% Inputs:
%    samples        - [NxD] - N samples with dimensionality D
%    c              - struct - user constraints
%
% Outputs:
%    drift               - user selection drift of samples
%    simspacePredictions - predicted similarity space locations of samples
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

drift = zeros(size(samples,1),1);

[classDistances,simspacePredictions] = applyConstraints(samples, c);
classBinary = false(size(classDistances,1),1);
classBinary(c.selectedClasses) = 1;
if sum(classBinary)==1
    distSELECTED = classDistances(classBinary,:);
else
    distSELECTED = min(classDistances(classBinary,:));
end
if sum(~classBinary)==1
    distNONSELECTED = classDistances(~classBinary,:);
else
    distNONSELECTED = min(classDistances(~classBinary,:));
end

if ~isempty(distNONSELECTED) % If not all classes are selected    
    drift = distSELECTED./(distSELECTED+distNONSELECTED); 
    drift = reshape(drift,length(drift),1); % Make sure it is a col vector    
end
end

