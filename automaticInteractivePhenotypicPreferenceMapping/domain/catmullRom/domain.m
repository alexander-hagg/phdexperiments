function d = domain(dof,varargin)
%domain - Get domain configuration for free form deformation of an n-polygon
%
% Syntax:  d = domain(dof)
%
% Inputs:
%    dof        - [1] - Set number of degrees of freedom
%
% Outputs:
%    d          - struct - Domain description struct. Please refer to the inline
%                          comments to understand which variables are expected
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2019; Last revision: 02-Apr-2019
%
%------------- BEGIN CODE --------------
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle')); % Random number stream
warning('off', 'MATLAB:MKDIR:DirectoryExists');warning('off', 'MATLAB:polyshape:repairedBySimplify');

%% Encoding
% Default number of degrees of freedom in representation
if strcmp(dof,'default'); d.dof = 16; else; d.dof = dof; end
name = 'symmetry'; if nargin > 1; name = varargin{1}; end
d.nGPUs = 1; if nargin > 2; d.nGPUs = varargin{2}; end


% FFD base shape (circle)
t = 0:2*pi/(d.dof/2):2*pi; t(end) = [];
x1 = 0.5*cos(t); y1 = 0.5*sin(t);
d.base = [x1,y1];
d.getPhenotype = @(genomes) getPhenotypeFFD(genomes,d.base);

% Set domain ranges
axialBoundAdaptation = 0.1; %0.01
radialBoundAdaptation = 0.25; %0.25
d.ranges(:,1) = [axialBoundAdaptation*ones(d.dof/2,1);-radialBoundAdaptation*pi*ones(d.dof/2,1)];
d.ranges(:,2) = [ ones(d.dof/2,1); radialBoundAdaptation*pi*ones(d.dof/2,1)];
d.validate        = 'validate'; % Validation function that is called when creating new solutions. Can contain any constraints.
disp(['Parameter Bounds: ' num2str(radialBoundAdaptation) ' / ' num2str(axialBoundAdaptation)]);

d.resolution = 64;
d.phenoDistMult = 30; % Increase distance between shapes when visualizing

%% Misc
d.flipMap = true;
d.projectdir = ['/home/' getenv('USER') '/archetype/'];
d.workdir = ['/scratch/ahagg2s/ARCHETYPE/workfolderCLUSTER/'];
d.simStartID = 1;

%% Categorization
d.featureLabels             = {'VAE 1','VAE 2'};

%% Fitness function
d.fitnessRange              = [0 1];
% pointSymmetry pointAsymmetry
% uiCircle uiLarge uiTriangle uiCircleTriangle uiCircleTriangleOblong
% user zeroFitness
FITNESSFUNCTION             = 'pointSymmetry';  rmpath(genpath('domain/catmullRom/fitnessFunctions')); addpath(genpath(['domain/catmullRom/fitnessFunctions/' FITNESSFUNCTION]));
d.selectPenalty = 'inverseDistance'; %'inverseDistance' 'relativeDistance' 'relativeDistanceOnlyPenalizeConstraintViolation'

%------------- END CODE --------------
