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
axialBoundAdaptation = 0.01;
radialBoundAdaptation = 0.25;
d.ranges(:,1) = [axialBoundAdaptation*ones(d.dof/2,1);-radialBoundAdaptation*pi*ones(d.dof/2,1)];
d.ranges(:,2) = [ ones(d.dof/2,1); radialBoundAdaptation*pi*ones(d.dof/2,1)];
d.validate        = 'validate'; % Validation function that is called when creating new solutions. Can contain any constraints.
disp(['Parameter Bounds: ' num2str(radialBoundAdaptation) ' / ' num2str(axialBoundAdaptation)]);

d.resolution = 64;

%% Misc
d.flipMap = true;
d.projectdir = ['/home/' getenv('USER') '/archetype/'];
d.workdir = ['/scratch/ahagg2s/ARCHETYPE/workfolderCLUSTER/'];
d.simStartID = 1;

%% Categorization
d.categorize = @(geno,pheno,p,d) categorize(pheno,d);

d.featureLabels             = {'area','perimeter','enstrophy','umax'};
d.featureMin                = [0.01       1        0.15       0.05]; %0.08 
d.featureMax                = [0.6        4        1.1        0.20]; % 0.18
    
if strcmp(name,'lettuce')
    d.selectedFeatures          = [1,3];
else
    d.selectedFeatures          = [1,2];
end

%% Fitness function
d.fitnessRange              = [0 1];
% Weight for acquisition function during surrogate-assisted QD (SAIL)
% More complex-to-model domains will have a higher variance and thus might need a lower weight value.
d.varCoef                   = 20;

d.fitfunPointSymm = @(geno) fitnessPointSymmetry(geno,d);
d.fitfunAreaCirc = @(geno) fitnessMaxAreaMinCircumference(geno,d);
d.fitfunLettuce = @(geno) fitnessLettuce(geno,d);
if strcmp(name,'lettuce') || strcmp(name,'lettuce_SAIL')
    d.fitfun = d.fitfunLettuce;
    disp(['Running Lettuce with ' int2str(d.nGPUs) ' GPUs']);
else
    d.fitfun = d.fitfunPointSymm;
end


%------------- END CODE --------------
