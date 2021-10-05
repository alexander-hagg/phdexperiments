function p = defaultParamSet(varargin)
% defaultParamSet - loads default parameters for QD algorithm 
% (here: MAP-Elites with grid archive)
%
% Syntax:  p = defaultParamSet(ncores)
%               ncores - number of parallel workers (one child per worker for efficiency)
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
p.nGens                     = 2^7;       % number of generations
p.nChildren                 = 2^4;      % number of children per generation
p.mutSigma                  = 0.1;      % mutation drawn from Gaussian distribution with this \sigma
p.resolution                = 16;  % Resolution of the map (in cells per dimension). 
p.extraMapValues            = {'uncertainty'}; % extra map values used in map struct
p.convergeLimit             = 0.05;
p.selectProcedure           = 'random'; %random or curiosity

% Visualization and management of saved data 
p.display.illu              = false;
p.display.illuMod           = 256;
p.numMaps2Save              = 100;

end

%------------- END CODE --------------
