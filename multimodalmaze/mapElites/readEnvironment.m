p.nIters                    = str2num(getenv('CFG_NITERS')); if isempty(p.nIters); p.nIters  = 2; end
p.startIter                 = str2num(getenv('CFG_STARTITER')); if isempty(p.startIter); p.startIter  = 2; end
p.selectionCriterionID      = str2num(getenv('CFG_SELCRITID')); if isempty(p.selectionCriterionID); p.selectionCriterionID = 2; end% Default: no selection
p.selectionValue            = str2num(getenv('CFG_SELVAL')); if isempty(p.selectionValue); p.selectionValue = 1; end
p.adjustObjective           = logical(str2num(getenv('CFG_ADJUSTOBJECTIVE'))); if isempty(p.adjustObjective); p.adjustObjective = true; end
p.loadFile                  = logical(str2num(getenv('CFG_LOADFILE'))); if isempty(p.loadFile); p.loadFile  = false; end
p.inputPath                 = getenv('CFG_INPUTPATH'); if isempty(p.inputPath); p.inputPath  = '/scratch/ahagg2s/GECCO2019/planner/baselines/'; end
p.inputFile                 = [char(getenv('CFG_INPUTFILE')) '.mat']; if length(p.inputFile)==4; p.inputFile  = '150.mat'; end
p.outputPath                = getenv('CFG_OUTPUTPATH'); if isempty(p.outputPath); p.outputPath  = '/scratch/ahagg2s/GECCO2019/runs/debugging/'; end
p.selectionMethod           = getenv('CFG_SELECTIONMETHOD'); if isempty(p.selectionMethod); p.selectionMethod='individual';end % 'class' or 'individual'
p.useDimReduction           = logical(str2num(getenv('CFG_DIMREDUCTION'))); if isempty(p.useDimReduction); p.useDimReduction=true;end
p.doSeeding                 = logical(str2num(getenv('CFG_SEEDING'))); if isempty(p.doSeeding); p.doSeeding=false;end

p.mutSigma                  = str2num(getenv('CFG_MUTRATE'));  if isempty(p.mutSigma); p.mutSigma=0.05;end
% Adjust mutation based on domain range
p.mutSigma = p.mutSigma*range(d.ranges);

d.penaltyWeight             = str2num(getenv('CFG_PWEIGHT')); if isempty(d.penaltyWeight); d.penaltyWeight = 0.5; end

disp(['Mutation set to: ' num2str(p.mutSigma) ' - or ' num2str(100*p.mutSigma/range(d.ranges)) '%']);
disp(['Penalty Weight set to: ' num2str(d.penaltyWeight)]);

