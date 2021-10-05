function [feature] = parsec_Categorize(samples, d)
%af_Categorize - Returns feature values between 0 and 1 for each dimension
%
% Syntax:  [feature] = af_Categorize(samples, d)
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
if d.alignedMap
    feature = samples(:,[2 3]);
else
    shape = d.express(samples);
    area = squeeze(polyarea(shape(1,:,:), shape(2,:,:)));
    camber = squeeze(shape(2,1:end/2,:)+shape(2,end:-1:end/2+1,:))/2;
    maxcamber = max(camber)';
    feature = [area maxcamber];
end

feature = (feature-d.featureMin)./(d.featureMax-d.featureMin);
feature(feature>1) = 1; feature(feature<0) = 0;

%------------- END OF CODE --------------
