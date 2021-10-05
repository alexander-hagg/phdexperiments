function d = domain(dof)
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
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle')); % Random number stream
warning('off', 'MATLAB:MKDIR:DirectoryExists');warning('off', 'MATLAB:polyshape:repairedBySimplify');

% Default number of degrees of freedom in representation
if strcmp(dof,'default')
    d.dof = 12;
else
    d.dof = dof;
end

% Features
d.categorize                = 'categorize';
d.featureLabels             = {'area','perimeter','maxspan','minspan','random'};
% Feature ranges
d.featureMin                = [0   1   0.5      0      0];
d.featureMax                = [3   7     3.5      0.05   1];
d.selectedFeatures          = [1    2]; % Default selection of features
d.nDims                     = 2; % Feature map resolution (do not change, other nD maps not supported as of yet)

d.debug                     = false;

d.baseShape = 'circle'; %'rectangle','circle'

% FFD base shape: CIRCLE.
t = 0:2*pi/(d.dof/2):2*pi;
t(end) = [];
x1 = 0.5*cos(t);
y1 = 0.5*sin(t);
d.baseCircle = [x1,y1];

% FFD base shape: SQUARE.
t = 0:2/(d.dof/2):2;
t(end) = [];
x1 = fx(t);
y1 = fy(t);
d.baseSquare = [x1,y1];

if strcmp(d.baseShape,'square')
    d.base = d.baseCircle;
elseif strcmp(d.baseShape,'circle')
    d.base = d.baseSquare;
end


% Fitness function
d.fitfun                    = @(geno) npolyObjective(geno,d);
d.fitnessRange              = [0 1];
% Weight for acquisition function during surrogate-assisted QD
d.varCoef                   = 20; % This might be domain dependent.
% More complex-to-model domains will have a higher variance and thus might need a lower weight value.

% Genotypic ranges
d.ranges          = [-0.5 0.5];
d.validate        = 'validate'; % Validation function that is called when creating new solutions. Can contain any constraints.

% Visualization
% This spacer is used when displaying multiple phenotypes in the same plot
d.spacer = 2.5;
d.flipMap = true;

%% OPTIONAL
% Domain explanation tab content, can contain description text that is shown to the user
d.description{1} = ['N-poly domain'];
d.description{2} = '';
d.description{3} = ['Description:      n-polygons, with n currently set to ' int2str(d.dof)];
d.description{4} = ['Representation:   sequence of (x,y) pairs.'];
d.description{5} = ['Fitness function: maximize radial symmetry.'];
% Show a high fitness and a low fitness phenotype
d.sampleHighFit         = zeros(1,d.dof);


function x = fx(t)
    x = zeros(1,length(t));
    x(t<=2.0) = -0.5;
    x(t<=1.5) = 0.5 - (t(t<=1.5)-1)/0.5;
    x(t<=1.0) = 0.5;
    x(t<=0.5) = -0.5 + t(t<=0.5)/0.5;


function y = fy(t)
    y = zeros(1,length(t));
    y(t<=2.0) = -0.5 + (t(t<=2.0)-1.5)/0.5;
    y(t<=1.5) = -0.5;
    y(t<=1.0) = 0.5 - (t(t<=1.0)-0.5)/0.5;
    y(t<=0.5) = 0.5;

%------------- END CODE --------------