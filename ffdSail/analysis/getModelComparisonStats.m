function [stats, p, d, collect] = getModelComparisonStats(experiments, p, d, names, varargin)
parse = inputParser;
parse.addRequired('experiments');
parse.addRequired('p');
parse.addRequired('d');
parse.addRequired('names');
parse.addOptional('percentiles',[25,75]);
parse.addOptional('bindiff',[0.05]);

parse.parse(experiments, p, d, names, varargin{:});
stats.cfg.percentiles           = parse.Results.percentiles;
stats.cfg.bindiff   = parse.Results.bindiff;

stats.cfg.xpNames = names;

pM = experiments{1}.predMap;
pM = ~cellfun(@isempty,{pM.cD});
[~,stats.cfg.iItr] = find(pM);
if isfield(stats.cfg, 'overrideShownPeriod')
    stats.cfg.iItr = stats.cfg.iItr(1:stats.cfg.overrideShownPeriod);
end
stats.cfg.nItr= size(stats.cfg.iItr,2);

mapsize(1) = size(experiments{1}.predMap(end).fitness,1);
mapsize(2) = size(experiments{1}.predMap(end).fitness,2);

stats.cfg.nTotalSamples = p.nTotalSamples;
[~, ~, segmentMap] = models_LoadCfg(d);
cfg.express = segmentMap{1};cfg.categorize = segmentMap{2};cfg.featureMin = segmentMap{3};cfg.featureMax = segmentMap{4};cfg.alignedMap = segmentMap{5};
stats.cfg.getCoordinates = @(x) feval(d.categorize, x, cfg);


% Maximum true bin fitness values in prediction map
stats.summary.maxFitMap = stats_maxFitMap(experiments);

% Find unique experiment configurations by name
[stats.cfg.xpUNIQ,~,stats.cfg.xpID] = unique(stats.cfg.xpNames);

% Collect data from experiment structs
for i=1:length(stats.cfg.xpID)
    for ii=1:length(experiments{i}.predMap)
        collect.raw.pred.fitness_true(i,ii,:)     = experiments{i}.predMap(ii).fitness_true(:);
        collect.raw.pred.fitness_true_rel(i,ii,:) = experiments{i}.predMap(ii).fitness_true(:)./stats.summary.maxFitMap(:);
        collect.raw.pred.invalid(i,ii,:)          = isnan(experiments{i}.predMap(ii).fitness_true(:));
    end
    for jj=1:length(experiments{i}.acqMap)
        collect.raw.acq.confidence(i,jj,:)        = experiments{i}.acqMap(jj).confidence(:);
    end
end

% Summarize experiments
for j=1:length(stats.cfg.xpUNIQ)
    collect.pred.bins.fitness_true(j,:,:)      = squeeze(nanmedian(collect.raw.pred.fitness_true(stats.cfg.xpID==j,:,:),1));
    collect.pred.bins.fitness_true_rel(j,:,:)  = squeeze(nanmedian(collect.raw.pred.fitness_true_rel(stats.cfg.xpID==j,:,:),1));
    collect.pred.bins.invalid(j,:,:)           = squeeze(nanmean(collect.raw.pred.invalid(stats.cfg.xpID==j,:,:),1));
    collect.acq.bins.confidence(j,:,:)         = squeeze(nanmedian(collect.raw.acq.confidence(stats.cfg.xpID==j,:,:),1));
    
    stats.pred.fitness_true.med(j,:)     = nanmedian(squeeze(collect.pred.bins.fitness_true(j,:,:)),2);
    stats.pred.fitness_true.prc(1,j,:)   = prctile(squeeze(collect.pred.bins.fitness_true(j,:,:)),stats.cfg.percentiles(1),2);
    stats.pred.fitness_true.prc(2,j,:)   = prctile(squeeze(collect.pred.bins.fitness_true(j,:,:)),stats.cfg.percentiles(2),2);
    
    stats.pred.fitness_true_rel.med(j,:) = nanmedian(squeeze(collect.pred.bins.fitness_true_rel(j,:,:)),2);
    stats.pred.fitness_true_rel.prc(1,j,:) = prctile(squeeze(collect.pred.bins.fitness_true_rel(j,:,:)),stats.cfg.percentiles(1),2);
    stats.pred.fitness_true_rel.prc(2,j,:) = prctile(squeeze(collect.pred.bins.fitness_true_rel(j,:,:)),stats.cfg.percentiles(2),2);
    
    stats.pred.invalid.mean(j,:)          = sum(squeeze(collect.pred.bins.invalid(j,:,:)),2);
    
    stats.acq.confidence.med(j,:)       = nanmedian(squeeze(collect.acq.bins.confidence(j,:,:)),2);
    stats.acq.confidence.prc(1,j,:)     = prctile(squeeze(collect.acq.bins.confidence(j,:,:)),stats.cfg.percentiles(1),2);
    stats.acq.confidence.prc(2,j,:)     = prctile(squeeze(collect.acq.bins.confidence(j,:,:)),stats.cfg.percentiles(2),2);
end

%% Compare experiments bin-wise
if length(stats.cfg.xpUNIQ) > 1
    for lo=1:length(stats.cfg.xpUNIQ)
        leaveout = lo;
        selmaps = setxor(1:length(stats.cfg.xpUNIQ),leaveout);
        for j = 1:length(selmaps)
            i = selmaps(j);
            %relativeFitness(i,lo,:) = (squeeze(collect.pred.bins.fitness_true(i,end,:)) - squeeze(collect.pred.bins.fitness_true(lo,end,:))) ./abs(squeeze(collect.pred.bins.fitness_true(lo,end,:)));
            relativeFitness(i,lo,:) = (squeeze(collect.pred.bins.fitness_true(i,end,:)) - squeeze(collect.pred.bins.fitness_true(lo,end,:))) ./stats.summary.maxFitMap(:);
            
            % Collect statistics
            stats.pred.bins.better(i,lo,:)   = relativeFitness(i,lo,:) > stats.cfg.bindiff;
            stats.pred.bins.worse(i,lo,:)    = relativeFitness(i,lo,:) < -stats.cfg.bindiff;
            stats.pred.bins.improved(i,lo,:) = stats.pred.bins.better(i,lo,:) - stats.pred.bins.worse(i,lo,:);
            stats.pred.bins.invalid(i,lo,:)  = isnan(relativeFitness(i,lo,:));
            
            stats.summary.bins.better(i,lo)   = sum(relativeFitness(i,lo,:) > stats.cfg.bindiff);
            stats.summary.bins.worse(i,lo)    = sum(relativeFitness(i,lo,:) < -stats.cfg.bindiff);
            stats.summary.bins.improved(i,lo) = -1*sum(relativeFitness(i,lo,:) < -stats.cfg.bindiff) + ...
                sum(relativeFitness(i,lo,:) > stats.cfg.bindiff);
            
            stats.pred.bins.relativeFitness(i,lo,:) = relativeFitness(i,lo,:);
        end
        stats.summary.bins.sum = sum(stats.summary.bins.improved,2);
    end
end
%% Get sample location densities
stats.samples.mapLocations = nan(length(stats.cfg.xpID),d.featureRes(1),d.featureRes(2));
for i=1:length(stats.cfg.xpID)
    stats.samples.genes(i,:,:) = squeeze(experiments{i}.samples.genes);
    stats.samples.mapLinIndx(i,:) = getBins( experiments{i}.samples.genes, stats.cfg.getCoordinates, d, experiments{i}.predMap(end).edges);
    for j=1:size(stats.samples.mapLinIndx,2)
        if isnan(stats.samples.mapLocations(i,stats.samples.mapLinIndx(i,j) )); stats.samples.mapLocations(i,stats.samples.mapLinIndx(i,j) ) = 0;end
        stats.samples.mapLocations(i,stats.samples.mapLinIndx(i,j) ) = stats.samples.mapLocations(i,stats.samples.mapLinIndx(i,j) ) + 1;
    end
end
%% Get mean map locations
for i=1:length(stats.cfg.xpUNIQ)
    stats.samples.meanMapLocations(i,:,:) = nanmean(stats.samples.mapLocations(stats.cfg.xpID==i,:,:));
end

%% Collect set of all optima and set of all optima+samples

disp(['Collect all optima']);
optimaLocations = [];
for i=1:length(experiments)
    loc = reshape(squeeze(experiments{i}.predMap(end).genes),mapsize(1)*mapsize(2),10);
    nanSolutions{i} = any(isnan(loc)');
    loc = loc(~any(isnan(loc)'),:);
    if i==1
        optimaIDs{i}= [1:size(loc,1)];
    else
        optimaIDs{i}= [optimaIDs{i-1}(end)+1:optimaIDs{i-1}(end)+size(loc,1)];
    end
    optimaLocations = [optimaLocations; loc];
end
%optimaLocationsAndSamples = optimaLocations;
%for i=1:length(experiments)
%    optimaLocationsAndSamples = [optimaLocationsAndSamples; experiments{i}.samples.genes];
%end
%nanSolutions{end+1} = logical(zeros(1,size(experiments{i}.samples.genes,1)*length(experiments)));

stats.agg.optimaLocations = optimaLocations;
stats.agg.optimaIDs = optimaIDs;
stats.agg.nanSolutions = nanSolutions;
%stats.agg.optimaLocationsAndSamples = optimaLocationsAndSamples;
stats.agg.nanSolutions = nanSolutions;

end