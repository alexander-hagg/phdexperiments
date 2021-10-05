function p = defaultParamSet(varargin)
% defaultParamSet - loads default parameters for MAP-Elites algorithm
%
% Syntax:  p = defaultParamSet(ncores)
%               ncores - number of parallel workers (one child per worker for efficiency)
%
% Outputs:
%   p      - struct - parameter struct
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018

%------------- BEGIN CODE --------------
p.nChildren                 = 30;
if nargin > 0; d = varargin{1}; end
% p.nChildren = d.ncores;
disp(['Number of children: ' int2str(p.nChildren)]);

p.mutSigma                  = 0.1;

p.constraintType            = 'prior'; %'none', 'prior' 'posterior'
p.nGens                     = 2^12;
p.numInitSamples            = 200;
p.selectProcedure           = 'random'; %random or curiousness

% Visualization and data management
p.display.illu              = false;
p.display.illuMod           = 25;
p.numMaps2Save              = 10;

% Read from environment variables
readEnvironment;

disp(['Selection method: ' p.selectionMethod]);

if p.useDimReduction
    disp('Using dimensionality reduction'); 
else
    disp('Not using dimensionality reduction'); 
end

p.selectionThreshold        = 0.99;
p.minSelectedSamples        = 50;
p.keepNonselectSeeds        = false;
p.dimreduxMethod			= 'tSNE';

end


