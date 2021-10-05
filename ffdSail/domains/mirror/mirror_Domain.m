function d = mirror_Domain(varargin)
%mirror_Domain - mirror Domain Parameters
% Returns struct with default for all settings of a car mirror domain
% including hyperparameters, and strings indicating functions for
% representation and evaluation. Direct parameter or free form deformation
% encodings can be chosen.
%
% Syntax:  d = mirror_Domain('encoding',ENCODING,'nCases',NCASES);
%
% Example:
%    d = mirror_Domain('encoding','ffd'  ,'nCases',10)
%    output = sail(sail,d);
%    d = mirror_Domain('encoding','param','nCases',2)
%    output = sail(sail,d);
%
%
% See also: sail, runSail

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: alexander.hagg@h-brs.de
% Dec 2017; Last revision: 11-Dec-2017
%------------- Input Parsing ------------
parse = inputParser;
parse.addOptional('encoding', 'ffd');
parse.addOptional('nCases'  , 1);
parse.addOptional('hpc', false);
parse.addOptional('lowres', false);
parse.addOptional('mirrorOnly', true);
parse.addOptional('runFolder', '/home/alex/Documents/phd/repositories/experiments/modules/prodigi');

parse.parse(varargin{:});
d.encoding    = parse.Results.encoding;
d.nCases      = parse.Results.nCases;
d.hpc         = parse.Results.hpc;
d.lowres      = parse.Results.lowres;
d.mirrorOnly  = parse.Results.mirrorOnly;
d.runFolder   = parse.Results.runFolder;

%%------------- BEGIN CODE --------------
disp(['Mirror domain: please ignore fitness values and use drag force values for analysis']);
d.namebase = 'mirror';
d.name = [d.namebase '_' d.encoding];
addpath(genpath(d.runFolder));
rmpath(genpath('ffdSail/domains'));
addpath(genpath('ffdSail/domains/mirror/'));

% - Scripts
% Common to any representations
d.preciseEvaluate   = 'mirror_PreciseEvaluate';    %
d.categorize        = 'mirror_Categorize';         %
d.createAcqFunction = 'mirror_CreateAcqFunc';      %
d.validate          = [d.name '_Validate'];

% - Genotype to Phenotype Expression
% Any representation should produce a fv struct and NX3 meshpoints
d.dof = 51;
base = 0.5+zeros(1,d.dof);
switch d.encoding
    case 'ffd'
        [~, ~, d.FfdP] =  mirror_ffd_Express(base,'ffdSail/domains/mirror/ffd/mirrorBase.stl');
        d.express  = @(x) mirror_ffd_Express(x, d.FfdP);
        [~, ~, d.FfdP_CFD] =  mirror_ffd_Express(base,'ffdSail/domains/mirror/ffd/mirrorBase_CFD.stl');
        d.express_CFD  = @(x) mirror_ffd_Express(x, d.FfdP_CFD);
        d.varCoef  = 1e0; % variance weight
    otherwise
        warning('No valid encoding defined');
end
d.base.fv   = d.express(base);
d.base.mesh = d.base.fv.vertices;

% - Alternative initialization [ TODO ]
d.loadInitialSamples = false;
d.initialSampleSource= '#notallFfdmirrors.mat';

% - Feature Space
% Map borders
d.featureRes = [16 16];
d.nDims      = length(d.featureRes);
d.featureLabels = {'TotalCurvature', 'RelativeLength', 'MirrorSurface'};
d.featureSelection = [1 2]; % Default selection (Total Curvature and Relative Length)
d.extraMapValues = {'dragForce','confidence'};

%% Features (domain specific)
%   Feature Borders
d.featureMin = [0.27  82  0];
d.featureMax = [0.35  120 40000];
d.featureMin = d.featureMin(d.featureSelection);
d.featureMax = d.featureMax(d.featureSelection);

% IDs of changeable vertices of base mesh
load('baseSubMeshIds.mat');d.features.subMeshIds = baseSubMeshIds;clear baseSubMeshIds;
% Vertice IDs of lines on which curvature is measured
load('verticeIDs.mat');d.features.curvature.ids = verticeIDs;clear verticeIDs;
% Reflective Surface
load('mirrorIDs.mat');d.features.mirror.ids = mirrorIDs;clear mirrorIDs;

%%

% - GP Models
modelCfgs = models_LoadCfg(d);
d.modelParsAcq(1,:) = modelCfgs(1,:);d.modelParsPre(1,:) = modelCfgs(1,:);

% Acquisition function
d.muCoef  = 1;                  % mean weight

% - Visualization
color8 = parula(8);
d.view = @(x) patch('Faces',x.faces,'Vertices',x.vertices,...
    'FaceColor',color8(3,:),'FaceAlpha',0.35);

%% Precise Evaluation
d.caseStart = 1;
d.nVals = 1;                    % # of values of interest, e.g. dragForce (1), or cD and cL (2)
d.maxDragForce = 25;            % Only use this to prevent really bad shapes to influence the surrogate model
d.minDragForce = 0.001;         % Only use this to prevent really bad shapes to influence the surrogate model


if d.hpc
    %% Cluster
    % % Cases are executed and stored here (cases are started elsewhere)
    d.openFoamFolder = ['/scratch/ahagg2s/sailCFD/'];
    if d.lowres
        d.openFoamTemplate = [d.runFolder '/ffdSail/domains/mirror/pe/ofTemplates/hpc1_short'];
    elseif d.mirrorOnly
        d.openFoamTemplate = [d.runFolder '/ffdSail/domains/mirror/pe/ofTemplates/hpc1_mirror_only'];        
    else
        d.openFoamTemplate = [d.runFolder '/ffdSail/domains/mirror/pe/ofTemplates/hpc1'];
    end
    % - There should be a folder called 'case1, case2, ..., caseN in this
    % folder, where N is the number of new samples added every iteration.
    % - Each folder has a shell script called 'caserunner.sh' which must be
    % run, this will cause an instance of openFoam to run whenever a signal is
    % given (to allow jobs to be run on nodes with minimal communication).
    
else
    %% Local
    % TODO: make script that creates folders and runs caserunners
    d.openFoamFolder = ['/scratch/ahagg2s/sailCFD/'];
    %if d.lowres
    %    d.openFoamTemplate = [d.runFolder '/ffdSail/domains/mirror/pe/ofTemplates/local_short'];
    %elseif d.mirrorOnly
        d.openFoamTemplate = [d.runFolder '/ffdSail/domains/mirror/pe/ofTemplates/local_mirror_only'];
    %else
    %    d.openFoamTemplate = [d.runFolder '/ffdSail/domains/mirror/pe/ofTemplates/local'];
    %end
end


%% Cases are executed and stored here
% d.caseRunner = [d.runFolder '/ffdSail/domains/mirror/pe/startCaseRunners.sh'];
% system(['cp ' d.caseRunner ' ' d.openFoamFolder]);
% for iCase = 1:d.nCases
%     system(['rm -rf ' d.openFoamFolder 'case' int2str(iCase)]);
%     system(['mkdir ' d.openFoamFolder 'case' int2str(iCase)]);
%     system(['cp -r ' d.openFoamTemplate '/* ' d.openFoamFolder 'case' int2str(iCase)]);
%     system(['cp ' d.caseRunner ' ' d.openFoamFolder ]);
% end

% %------------- END OF CODE --------------






