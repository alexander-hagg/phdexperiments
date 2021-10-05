function d = domain_Maze(varargin)
%DOMAIN_MAZE Get domain configuration for multimodal maze
%   Radius of maze = 336 pixels
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle')); % Random number stream
warning('off', 'MATLAB:MKDIR:DirectoryExists');

d.nDims                     = 2;
d.mazeFileName              = 'mediumRound';
d.featureMin                = 0;
d.featureMax                = 400;
d.categorize                = 'categorize';
d.extraMapValues            = {'alignment','ring1','ring2','ring3','orgFitness','penalty'};
d.featureLabels             = {'x','y'};
d.featureRes                = [30 30];
d.center                    = [200 200];
d.startPosition             = d.center;
d.ringDiameter(1)           = 103;
d.ringDiameter(2)           = 192;
d.ringDiameter(3)           = 277;
d.goalRings{1}              = [200,150;243,228;161,225];
d.goalRings{2}              = [282,150;202,296;120,150];
d.goalRings{3}              = [200,60;320,275;85,275];
d.debug                     = false;

d.penaltyWeight             = str2num(getenv('CFG_PWEIGHT')); if isempty(d.penaltyWeight); d.penaltyWeight = 0.5; end
disp(['Penalty weight: ' num2str(d.penaltyWeight)]);

d.tmpdir = getenv('JOBTMPDIR'); if isempty(d.tmpdir); d.tmpdir='/tmp';end
disp(['tmp dir: ' d.tmpdir]);
mkdir(d.tmpdir);

%% Reconfigure maze
% d.ncores = str2num(getenv('CFG_NCORES'))-2; if isempty(d.ncores); d.ncores      = 3; end
% if d.ncores > 3; pc = parcluster('local'); pc.JobStorageLocation = strcat(getenv('JOBTMPDIR')); parpool(pc, d.ncores);end
d.theta       = str2num(getenv('CFG_THETA')); if isempty(d.theta); d.theta      = 30; end
fileExt       = string(d.theta);
d.mazeCfgBaseFile = d.mazeFileName;
xmlOutFile = simSetAngle(d.mazeCfgBaseFile, d.theta, fileExt);
d.mazeCfgFile = xmlOutFile;

d.representation = getenv('CFG_REPRESENTATION'); if isempty(d.representation); d.representation = 'controller';end %'planner' 'controller'
disp(['Using ' d.representation ' domain']);

if strcmp(d.representation,'planner')
    %% Planner
    d.ranges      = [-200 200];
    d.dof         = 14;
    d.phenotypeLength = 2000; % Maximal phenotype length
    d.validate    = 'validateChildren_plan';
    d.map = importdata(['domains/mazeDomain/simulator/worlds/' d.mazeCfgBaseFile '.pbm']);
    d.evalFcn     = @(samples) eval_maze_plan(samples, d, false);
    d.startPosition = d.startPosition + 10*(rand(size(d.startPosition))-0.5); % Randomize starting position a bit
    d.flipMap = false;
    
elseif strcmp(d.representation,'controller')
    
    %% Simulator/ANN
    d.useGoalQuadrantSensors    = true;
    d.useRNN      = logical(str2num(getenv('CFG_RNN'))); if isempty(d.useRNN); d.useRNN = false; end
    if d.useRNN
	disp('Using RNN');
    else
	disp('Not using RNN');
    end
    d.numHidden   = str2num(getenv('CFG_NUMHIDDEN')); if isempty(d.numHidden); d.numHidden  = 5; end
    
    d.ncores   = str2num(getenv('CFG_NCORES')); if isempty(d.ncores); d.ncores  = 4; end
    disp(['Number of CPU cores available: ' int2str(d.ncores)]);
    
    if d.useGoalQuadrantSensors
        disp('Running with goal quadrant sensors'); d.numInputs = 8;d.useGoalQuadrantSensors= true;
    else
        disp('Running without goal quadrant sensors');d.numInputs = 4; d.useGoalQuadrantSensors= false;
    end
    
    if ~d.useRNN
        disp('Feed forward controllers (MLPs)');
        d.dof                       = (d.numInputs) * d.numHidden + (d.numHidden + 1) * 2;
    else
        disp('Recurrent controllers (Elman networks)');
        d.dof                       = (d.numInputs) * d.numHidden + ...   % input to hidden (full)
            d.numHidden * 2 + ...             % hidden to output (full)
            d.numHidden + ...                 % hidden to context (1-1)
            d.numHidden + ...                 % context to itself (1-1)
            d.numHidden * d.numHidden + ...   % context to hidden (full)
            d.numHidden + ...                 % bias context
            2;                                % bias outputs
    end
    
    d.timesteps       = 1000;
    d.ranges          = [-3 3];
    d.evalFcn         = @(samples) eval_maze(samples, d, false);
    d.phenotypeLength = d.timesteps;
    d.validate        = 'validateChildren';
    d.flipMap = true;    

end

%% Individual's genome and phenotype
d.sampleInd.genome    = nan(d.dof,1);
d.sampleInd.phenotype = nan(3*d.phenotypeLength,1);

%% Objective Function: Path length
d.metricFitness               = @metricFitness;
d.fitfun                      = @(X) objective(X, d.evalFcn, d.metricFitness, [], d.penaltyWeight);

end

