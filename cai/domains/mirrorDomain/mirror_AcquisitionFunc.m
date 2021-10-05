function [fitness,predValue] = mirror_AcquisitionFunc(drag,d)
%mirror_AcquisitionFunc - Infill criteria based on uncertainty and fitness
%
% Syntax:  [fitness, dragForce] = mirror_AcquisitionFunc(drag,d)
%
% Inputs:
%   drag -    [2XN]    - dragForce mean and variance
%   d    -             - domain struct 
%   .varCoef  [1X1]    - uncertainty weighting for UCB
%
% Outputs:
%    fitness   - [1XN] - Fitness value (lower is better)
%    predValue - [1XN] - Predicted drag force (mean)
%

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 01-Aug-2017

%------------- BEGIN CODE --------------

fitness = (drag(:,1)*d.muCoef) - (drag(:,2)*d.varCoef); % better fitness is lower fitness  
predValue{1} = drag(:,1);
predValue{2} = drag(:,2);

%------------- END OF CODE --------------
