function d = cppn_Domain
%cppn_Domain - CPPN-derived feed forward deformation Airfoil Domain Parameters 
%
%Returns struct with default for all settings of the FFD-airfoil domain
%including hyperparameters, and strings indicating functions for
%representation and evaluation.
%
% Syntax:  d = cppn_Domain;
%
% Example: 
%    output = sail(sail,cppn_Domain);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: sail, runSail

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: adam.gaier@h-brs.de
% Jun 2017; Last revision: 20-Aug-2017

%------------- BEGIN CODE --------------
d.name = 'CPPN';
rmpath( genpath('domains'));
addpath(genpath('domains/foilFFD/'));

% Scripts
d.preciseEvaluate   = 'ffd_PreciseEvaluate';
d.categorize        = 'ffd_Categorize';
d.createAcqFunction = 'ffd_CreateAcqFunc';
d.validate          = 'ffd_ValidateChildren';
d.saveData          = 'ffd_RecordData';

% Alternative initialization
d.loadInitialSamples = false;
d.initialSampleSource= '';

% Genotype to Phenotype Expression
d.dof       = 10;
d.express   = @(x) cppnRaeY(x);
d.base      = loadBaseAirfoil(d.express, d.dof);

% Feature Space
d.featureRes    = [25 25];
d.nDims         = length(d.featureRes);
d.featureMin    = [0.0440 0.1500];
d.featureMax    = [0.1588 0.5175];
d.featureLabels = {'Z_{up}','X_{up}'}; % {X label, Y label}

% GP Models
d.gpParams(1)= paramsGP(d.dof); % Drag
d.gpParams(2)= paramsGP(d.dof); % Lift

% Acquisition function
d.varCoef = 20; % variance weight
d.muCoef  = 1; % mean weight 

% Domain Specific Map Values
d.extraMapValues = {'cD','cL','confidence'};

%------------- END OF CODE --------------




















