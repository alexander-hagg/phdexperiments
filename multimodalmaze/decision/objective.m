function [adjustedFitness, trajectories, values] = objective(X, evalFcn, metricFitness, constraintSet, penaltyWeight, varargin)
%OBJECTIVE Generic objective function
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018
%
%------------- BEGIN CODE --------------  

if nargin > 5
    fitness = varargin{1};
    values = varargin{2};
else
    [trajectories,values] = evalFcn(X);
    fitness = metricFitness(trajectories);
end
        
penalty = zeros(size(X,1),1);
if ~isempty(constraintSet)
    for iT=1:length(constraintSet)
        if ~strcmp(constraintSet(iT).constraints,'none') && ~isempty(constraintSet(iT).constraints)
            penalty = penalty + constraintPenalty(X,constraintSet(iT).constraints);
        end
    end
end

values{end+1} = fitness';
values{end+1} = penalty';
adjustedFitness = fitness .* (1 + penalty.*penaltyWeight);
end

