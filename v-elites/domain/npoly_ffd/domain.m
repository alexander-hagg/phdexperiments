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
d.featureMin                = [0   1.75   0.5      0      0];
d.featureMax                = [4   15     3.5      0.05   1];
d.selectedFeatures          = [1    2]; % Default selection of features
d.nDims                     = 2; % Feature map resolution (do not change, other nD maps not supported as of yet)

d.debug                     = false;

% FFD base shape. Specific to this domain, as encoding uses a free form deformation of a base shape.
t = 0:2*pi/(d.dof/2):2*pi;
t(end) = [];
x1 = 0.5*cos(t);
y1 = 0.5*sin(t);
[theta,rho] = cart2pol(x1,y1);
d.base = [theta;rho];

% Fitness function
d.fitfun                    = @(X) npolyObjective(X,d);
d.fitnessRange              = [0 1];

% Genotypic ranges
d.ranges          = [-1 1];
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

end

%------------- END CODE --------------