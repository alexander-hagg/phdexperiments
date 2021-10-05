function output = produqd(p,d)
%PRODUQD Summary of this function goes here

% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Jan 2018; Last revision: 25-Jan-2018

%------------- BEGIN CODE --------------

if nargin==0; output = defaultParamSetProduqd; return; end

%% 0 - Initialize samples
if ~d.loadInitialSamples
    [observation, value] = initialSampling(d,p.qd.nInitialSamples);
    d.initialSampleSource = ['testSamples_' int2str(randi(100000)) '.mat'];
    save(d.initialSampleSource, 'observation', 'value');
    d.loadInitialSamples = true;
elseif ~isfield(p, 'sailInput')
    load(d.initialSampleSource); % Contains fields: 'observation', 'value'
    randomPick = randperm(size(observation,1),p.qd.nInitialSamples); %#ok<NODEF>
    observation = observation(randomPick,:);
    value = value(randomPick,:); %#ok<NODEF>
end

%% 1 - Design Iteration
if isfield(p,'inputRuns')
    output = p.inputRuns;
    startIter = length(output)+1;
else
    startIter = 1;
end

for it = startIter:p.numIterations
    
    %% 1a - Illuminate
    % Use SAIL to illuminate the feature map according to set parameter ranges
    % Check for SAIL run prior
    if it==1 && isfield(p, 'sailInput')
        disp(['Iteration ' int2str(it) ' | Running QD']);
        output{it}.sail = p.sailInput;
    else
        disp(['Iteration ' int2str(it) ' | Running QD']);
        output{it}.sail = sail(p.qd,d);
    end
    %% 1b - Get concepts and prototypes
    if it==1 && isfield(p, 'simspaceInput')
        output{it}.optima = p.simspaceInput.optima;
        output{it}.conceptLabels = p.simspaceInput.conceptLabels;
        output{it}.latent = p.simspaceInput.latent;
        output{it}.concepts = p.simspaceInput.concepts;
        output{it}.prototypes = p.simspaceInput.prototypes;
        output{it}.prototypesLatent = p.simspaceInput.prototypesLatent;
    else
        disp(['Iteration ' int2str(it) ' | Extract Concepts']);
        output{it}.optima = reshape(output{it}.sail.predMap(end).genes,prod(p.qd.data.predMapRes),d.dof);
        output{it}.optima(any(isnan(output{it}.optima')),:) = [];
        
        % Extract concepts
        [output{it}.conceptLabels, output{it}.latent, ~, ~, ~, output{it}.concepts] = ...
            dimReducedClustering( output{it}.optima, 'tSNE', 2);
        
        % Extract prototypes, excluding non-assigned concept ID '0'
        for ii=1:max(output{it}.concepts.uniqid)
            samplesLatent = output{it}.latent(output{it}.conceptLabels==ii,:);
            samples = output{it}.optima(output{it}.conceptLabels==ii,:);
            [~,output{it}.prototypesLatent(ii,:),~,~,protIDs] = kmedoids(samplesLatent,1);
            output{it}.prototypes(ii,:) = samples(protIDs,:);
        end
    end
    
    
    %% 1c - Select concept
    disp(['Iteration ' int2str(it) ' | Automated Selection and Constraint Extraction']);
    if it > 1 && p.selectCriterion.oneTime
        criterion.type = 'nearest';criterion.valType = 'rank';criterion.value = [1];
        p.oldconceptID = output{1}.conceptSelection.id;
        p.qd.concept.id = selectConcept(output{it}.concepts.sizes(2:end), criterion, p.oldconceptID, output{1}.prototypesLatent, output{it}.prototypesLatent);
    else
        p.qd.concept.id = selectConcept(output{it}.concepts.sizes(2:end), p.selectCriterion, output{it}.prototypes);
        p.oldconceptID = p.qd.concept.id;
    end
    p.qd.concept.members = output{it}.optima(ismember(output{it}.conceptLabels,p.qd.concept.id),:);
    output{it}.conceptSelection = p.qd.concept;
    
    %% 1d - Reconfigure illumination
    % Train constraints if applicable
    if isfield(p.qd, 'constraints')
        if ~isnan(p.qd.constraints.threshold)
            if ~p.selectCriterion.oneTime || it == 1
                disp(['Training constraint model']);
                p.qd.constraintModel = trainConstraints(output{it}.optima,output{it}.latent);
                p.qd.oldconceptID = p.oldconceptID;
                p.qd.concept.allLabels = output{it}.conceptLabels;
            end
        end
    end
    
    % Save training samples
    d.initialSampleSource = ['PRODIGI_samples_' int2str(randi(100000)) '.mat'];
    observation = [output{it}.sail.tModelsAcq{end,1}.trainInput];
    value = nan(size(observation,1),size(output{it}.sail.tModelsAcq,2));
    for i = 1:size(output{it}.sail.tModelsAcq,2)
        value(:,i) = output{it}.sail.tModelsAcq{end,i}.trainOutput;
    end
    save(d.initialSampleSource,'observation','value');
    d.loadInitialSamples = true;
    nAddedSamples = p.qd.nTotalSamples - p.qd.nInitialSamples;
    p.qd.nInitialSamples    = size(observation,1); % Observation set might be smaller than expected, so just assign its length
    p.qd.nTotalSamples      = p.qd.nInitialSamples + nAddedSamples;
    
    output{it}.p            = p; % Save configuration
end


end

