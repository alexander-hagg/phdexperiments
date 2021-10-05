function [predMap, allMaps] = createPredictionMap(model,p,d, varargin)
%createPredictionMap - Produce prediction map from surrogate model
%
% Syntax:  predictionMap = createPredictionMap(gpModels,p,d)
%
% Inputs:
%    gpModels   - GP model produced by SAIL
%    p          - SAIL hyperparameter struct
%    d          - Domain definition struct
%
% Outputs:
%    predMap - prediction map
%    .fitness     [Rows X Columns]
%    .genes       [Rows X Columns X GenomeLength]
%    .'otherVals' [Rows X Columns]
%    percImproved    - percentage of children which improved on elites
%
% Example: 
%    p = sail;
%    d = af_Domain;
%    output = sail(d,p);
%    predMap = createPredictionMap(output.model,p,d,'featureRes',[50 50]);
%    viewMap(predMap.fitness,d, predMap.edges)
%
% Other m-files required: mapElites  nicheCompete updateMap d.createAcqFunction
%
% See also: sail, mapElites, runSail

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: adam.gaier@h-brs.de
% Jun 2017; Last revision: 03-Aug-2017
%
%------------- BEGIN CODE --------------

if nargin > 4; observation = varargin{1};disp('Initial samples from user');end
if nargin > 5; figHandleMap = varargin{2};end

% Construct functions
d.varCoef = 0; % no award for uncertainty
acqFunction = @(x) ucb(x, model, d, p);

% Seed map with precisely evaluated solutions
if ~exist('observation','var')
    disp('Initial samples from model');
    observation = model.trainInput;
end
    

[fitness,features] = acqFunction(observation);

predMap                                              = createMap(d, p);
[replaced, replacement, features, percImprovement]   = nicheCompete(observation, fitness, predMap, d, p, features);
predMap                                              = updateMap(replaced,replacement,predMap,fitness,observation, features);                    
    
    
% Illuminate based on surrogate models
[predMap,~,~,allMaps] = illuminate(predMap,acqFunction,p,d);
%------------- END OF CODE --------------