function p = defaultParamSetProduqd()
%DEFAULTPARAMSETPRODUQD Load default parameters for iProD algorithm
%
% Syntax:  p = defaultParamSetPRODIGI
%
% Outputs:
%   p      - struct - parameter struct
%

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2018; Last revision: 25-Jan-2018

%------------- BEGIN CODE --------------

p.qd                        = defaultParamSet;  % Load default Quality-Diversity hyperparameters

% Iterations
p.numIterations             = 2;                % Number of engineering design iterationss

% Concept Selection Types
p.selectCriterion.type      = 'size';           % 'nearest', 'size', variance', 'sparsity', 'map', 'prototypedist'
p.selectCriterion.valType   = 'rank';           % 'rank', 'percentile'
p.selectCriterion.value     = [1];              % [1 2 3] ('type' top three), [0.95 1.0] ('type' top 5%)
                                                % [-1 -2 -3] ('type' bottom three), [0.0 0.05] ('type' bottom 5%)
p.selectCriterion.oneTime   = true;            % If true, selection is only applied once, 
                                                % after which selection is based on proximity to prototype in first selection

%------------- END OF CODE --------------


end

