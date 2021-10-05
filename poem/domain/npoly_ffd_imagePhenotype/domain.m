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
% Nov 2019; Last revision: 13-Nov-2019
%
%------------- BEGIN CODE --------------
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle')); % Random number stream
warning('off', 'MATLAB:MKDIR:DirectoryExists');warning('off', 'MATLAB:polyshape:repairedBySimplify');

% Default number of degrees of freedom in representation
if strcmp(dof,'default')
    d.dof = 16;
else
    d.dof = dof;
end

% FFD base shape (circle)
    t = 0:2*pi/(d.dof/2):2*pi;
    t(end) = [];
    x1 = 0.5*cos(t);
    y1 = 0.5*sin(t);
d.base = [x1,y1];

d.getPhenotype = @(genomes) getPhenotypeFFD(genomes,d.base);


% Fitness function
d.fitfunPointSymm                    = @(geno) fitnessPointSymmetry(geno,d);
d.fitfunAreaCirc                    = @(geno) fitnessMaxAreaMinCircumference(geno,d);

d.fitfun = d.fitfunPointSymm;

d.fitnessRange              = [0 1];
% Weight for acquisition function during surrogate-assisted QD
d.varCoef                   = 20; % This might be domain dependent.
% More complex-to-model domains will have a higher variance and thus might need a lower weight value.

% Genotypic ranges
d.ranges(:,1) = [0*ones(d.dof/2,1);-pi*ones(d.dof/2,1)];
d.ranges(:,2) = [2*ones(d.dof/2,1);pi*ones(d.dof/2,1)];
d.validate        = 'validate'; % Validation function that is called when creating new solutions. Can contain any constraints.

% Visualization
d.flipMap = true;

%------------- END CODE --------------