function run_cai()
%RUN_PRODUQD Prototype Discovery Using Quality-Diversity
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 03-Jan-2019
%
%------------- BEGIN CODE --------------
%% Configure experiment
addpath(genpath('.'));
d = domain_Maze; p = defaultParamSet(d);

if p.loadFile
    disp(['Preloading MAP-Elites data from: ' p.inputPath p.inputFile]);
    loadedData = load([p.inputPath p.inputFile]); 
    runData = loadedData.runData; samples = loadedData.samples; fitness = loadedData.fitness; phenotype = loadedData.phenotype; values = loadedData.values;
    acqMap = loadedData.runData{1}.acqMap;percImproved = loadedData.runData{1}.percImproved;percValid = loadedData.runData{1}.percValid;percFilled = loadedData.runData{1}.percFilled;
    d = runData{1}.d; p = runData{1}.p; readEnvironment;
    
    if strcmp(d.representation,'controller')
        p.iter2.nGens = 2^12;
    elseif strcmp(d.representation,'planner')    
        p.iter2.nGens = 2^12;
    end
    p.iter2.selectProcedure = 'random';
    
    % Reset temporary directy
    d.tmpdir = getenv('JOBTMPDIR'); if isempty(d.tmpdir); d.tmpdir='/tmp';end
    disp(['tmp dir: ' d.tmpdir]);
    mkdir(d.tmpdir);
    if strcmp(d.representation,'planner')
        d.evalFcn     = @(samples) eval_maze_plan(samples, d, false);
    elseif strcmp(d.representation,'controller')
        d.evalFcn     = @(samples) eval_maze(samples, d, false);
    end
    d.fitfun = @(X) objective(X, d.evalFcn, d.metricFitness, [], d.penaltyWeight);
else
    if strcmp(d.representation,'controller')
        sobSequence                         = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');
        samples                             = 0.1*d.ranges(2)*((2*sobSequence(1:(1+p.numInitSamples)-1,:))-1);
    elseif strcmp(d.representation,'planner')
        samples = 0.1*randn(p.numInitSamples*10,d.dof)*d.ranges(2);
        samples(:,1:end/2) = samples(:,1:end/2) + (d.startPosition(1)-d.center(1));
        samples(:,end/2+1:end) = samples(:,end/2+1:end) + (d.startPosition(2)-d.center(2));
    end
    [fitness,phenotype,values]          = d.fitfun(samples);
    validInds                           = feval(d.validate, samples, d);
    disp(['Found ' int2str(sum(validInds)) '/' int2str(length(validInds)) ' valid initial samples']);
    samples = samples(validInds,:);
    fitness = fitness(validInds,:);
    phenotype = phenotype(validInds,:,:);
    for valI = 1:length(values); values{valI} = values{valI}(validInds);end
end

%%
for iT = p.startIter:(p.startIter+p.nIters-1)
    
    if ~(p.loadFile && iT == p.startIter)
        disp('MAP-Elites started');
        obsMap                              = createMap(d.featureRes, d.sampleInd, d.extraMapValues);
        [replaced, replacement]             = nicheCompete(samples, fitness, phenotype, obsMap, d);
        obsMap                              = updateMap(replaced,replacement,obsMap,fitness,samples,values,d.extraMapValues);
        
        %profile on
        [acqMap, percImproved, percValid, ~, allMaps, percFilled] = mapElites(d.fitfun,obsMap,p,d);
        %profile off
        runData_EXTRA{iT}.allMaps = allMaps; runData_EXTRA{iT}.d = d; runData_EXTRA{iT}.p = p;
        save([p.outputPath int2str(d.theta) '_allMaps.mat'], 'runData_EXTRA', '-v7.3'); disp('SAVED');
        
        samples = reshape(acqMap.genes,d.featureRes(1)*d.featureRes(2),d.dof); samples = samples(~any(isnan(samples)'),:); % Filter out NaNs
        
        runData{iT}.d = d;runData{iT}.p = p; runData{iT}.acqMap = acqMap; runData{iT}.percImproved = percImproved;runData{iT}.percValid = percValid;runData{iT}.percFilled = percFilled;
        disp('MAP-Elites done');
    end
    
    % TODO: get rid of this hack. Recalculate values...
    samples = samples(~any(isnan(samples)'),:); % Filter out NaNs
    [fitness, phenotype, values] = d.fitfun(samples);
    
    %% Extract classes
    disp('Classification started'); classification = extractClasses(samples,p.dimreduxMethod); disp('Classification and prototyping done');
    runData{iT}.classification = classification;
    
    if p.selectionCriterionID < 9999 && p.selectionValue < 9999 && (iT < (p.startIter+p.nIters-1))
        %% Train Selection
        [isSelected,runData{iT}] = selectionProcedure(values, classification, p, runData{iT});
        p.constraintSet(iT).constraints = runData{iT}.constraints;
        
        %% Set and recalculate, with or without seeding
        if p.doSeeding
            %% Set seeds
            disp('Setting seeds');
            d.fitfun = @(X) objective(X, d.evalFcn, d.metricFitness, [], d.penaltyWeight);
            fitness = fitness(isSelected); phenotype = phenotype(isSelected,:,:); for valI = 1:length(values); values{valI} = values{valI}(isSelected);end
            samples = samples(isSelected,:);
        end
        if p.adjustObjective
            %% Adjust objective function with user decision hypersurface model
            disp('Adjusting objective function to include user decision hypersurface model');
            d.fitfun = @(X) objective(X, d.evalFcn, d.metricFitness, p.constraintSet, d.penaltyWeight);
            [fitness, phenotype, values] = d.fitfun(samples);
            % Filter out NaNs
            fitness = fitness(~any(isnan(samples)')); phenotype = phenotype(~any(isnan(samples)'),:,:); for valI = 1:length(values); values{valI} = values{valI}(~any(isnan(samples)'));end
            samples = samples(~any(isnan(samples)'),:);
        end
        
    end
    
    save([p.outputPath int2str(d.theta) '.mat'], 'runData', 'samples', 'fitness', 'phenotype', 'values', 'd','p','-v7.3');
    
    p.nGens = p.iter2.nGens;
    p.selectProcedure = p.iter2.selectProcedure;
    
end

end
%------------- END CODE --------------
