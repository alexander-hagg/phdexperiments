function p = defaultParamSet(varargin)
% defaultParamSet - loads default parameters for QD algorithm 
% (here: MAP-Elites with Voronoi archive)
%
% Syntax:  p = defaultParamSet(varargin)
%
% Outputs:
%   p      - struct - QD configuration struct
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

p.numInitSamples            = 2^6;      % number of initial samples
p.nGens                     = 2^8;       % number of generations
p.nChildren                 = 2^5;      % number of children per generation
p.mutSigma                  = 0.2;      % mutation drawn from Gaussian distribution with this \sigma
%p.maxBins                   = 2^3;     % not used yet but intention is to
                                        % enable adjusting feature resolution 
                                        % based on number of requested bins to
                                        % control complexity
p.featureResolution         = 0.04;     % resolution of the map is controlled by local competition distance threshold

% Selection
p.penaltyWeight             = 2;        % User selection drift, weight for soft constraints
p.driftThreshold            = 0.5;      % User selection drift, threshold for hard user constraint

% Visualization and data management
p.display.illu              = true;
p.display.illuMod           = 25;
p.extraMapValues            = {'fitnessAdjustment','drift'};

end


%------------- END OF CODE --------------