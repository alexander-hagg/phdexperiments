function p = defaultParamSet(d)
% defaultParamSet - loads default parameters for SAIL algorithm
%
% Syntax:  p = defaultParamSet
%
% Outputs:
%   p      - struct - parameter struct
%

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 02-Aug-2017

%------------- BEGIN CODE --------------

% MAP-Elites Parameters
p.nChildren         = 32; 
p.mutSigma          = 0.1; 
p.nGens             = 512;

% Infill Parameters
p.nInitialSamples   = 8;
p.nAdditionalSamples= 8;
p.nTotalSamples     = 32;
p.trainingMod       = 3;

% Display Parameters
p.display.figs      = true;
p.display.gifs      = false;
p.display.illu      = true;
p.display.illuMod   = 100;%p.nGens;

% Data Gathering Parameters
p.data.outSave      = true;
p.data.outMod       = 10;
p.data.mapEval      = false;
p.data.mapEvalMod   = p.nTotalSamples;
p.data.mapEvalSteps = p.nInitialSamples:p.data.mapEvalMod:p.nTotalSamples;
p.data.modelSaveSteps   = p.data.mapEvalSteps;
p.data.outPath      = '';

p.retryInvalid      = true;

p.fitfun            = @(X) objective(X, d.evalFcn, [], d.penaltyWeight);
p.constraints       = [];

%------------- END OF CODE --------------

