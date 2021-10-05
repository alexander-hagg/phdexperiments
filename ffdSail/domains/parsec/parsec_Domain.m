function d = parsec_Domain(tmpdir)
%af_Domain - Airfoil Domain Parameters
%Returns struct with default for all settings of airfoil domain including
%hyperparameters, and strings indicating functions for representation and 
%evaluation.
%
% Syntax:  d = parsec_Domain;
%
% Example: 
%    output = sail(sail,af_Domain);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required:
%   /domains/airfoil/airfoilTools/demo/raeParsec_aft.mat
%   /domains/airfoil/airfoilTools/demo/raeRange_aft.mat
%
% See also: sail, runSail

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: adam.gaier@h-brs.de
% Jun 2017; Last revision: 01-Aug-2017

%------------- BEGIN CODE --------------
d.name = 'PARSEC';
rmpath( genpath('domains'));
addpath(genpath('domains/parsec/'));
d.tmpdir = tmpdir;

% Scripts
d.preciseEvaluate   = 'parsec_PreciseEvaluate';
d.categorize        = 'parsec_Categorize';
d.createAcqFunction = 'parsec_CreateAcqFunc';
d.validate          = 'parsec_ValidateChildren';

% Alternative initialization
d.loadInitialSamples = false;
d.initialSampleSource= '';

% Genotype to Phenotype Expression
d.dof = 10;
load('raeParsec_aft.mat'); load('raeRange_aft.mat');  
d.base = loadBaseAirfoil(raeParsec, raeRange, d.tmpdir);
d.express = setExpression(d.base.parsec, d.base.range, d.dof);
d.base.genome = [0.115455731631420,0.769456784342656,0.139077011276969,0.278798623701840,0.124359615964012,0.151578562448986,0.751881051352209,0.307595979711475,0.262943052636945,0.762983847665074];
d.fitnessPenalty_Area = false;

% Feature Space
d.featureRes = [25 25];
d.nDims      = length(d.featureRes);
d.alignedMap = false;
if d.alignedMap
    d.featureMin = [0 0];
    d.featureMax = [1 1];
    d.featureLabels = {'X_{up}','Z_{up}'};
else
    d.featureMin = [0.05 0];
    d.featureMax = [0.2 0.08];
    d.featureLabels = {'sum area','max. camber'};
end

% Acquisition function
d.varCoef = 20; % variance weight
d.muCoef  = 1; % mean weight 
d.prevMaxVar = []; % keep track of previous maximum variances

% Domain Specific Map Values
d.extraMapValues = {'cD','cL','confidence'};

%------------- END OF CODE --------------




















