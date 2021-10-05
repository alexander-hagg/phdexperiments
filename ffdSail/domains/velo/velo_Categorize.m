function [feature] = velo_Categorize(samples, d)
%velo_Categorize - Returns feature values between 0 and 1 for each dimension
%
% Syntax:  [feature] = velo_Categorize(samples, d)
%
% Inputs:
%   samples  - [N X genomeLength] - uncategorized solutions
%   d        -  Domain description struct
%   .featureMin [1X1]   - minimum feature value
%   .featureMax [1X1]   - maximum feature value
%
% Outputs:
%    feature - [MXN] - Feature values for each individual
%

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 02-Aug-2017

%------------- BEGIN CODE --------------

parfor i=1:size(samples,1)
    [~, ~, feature(i,:)] = expressVelo(samples(i,:),'categoryOnly',true);
end
feature = (feature-d.featureMin)./(d.featureMax-d.featureMin);
feature(feature>1) = 1; feature(feature<0) = 0;
%------------- END OF CODE --------------
