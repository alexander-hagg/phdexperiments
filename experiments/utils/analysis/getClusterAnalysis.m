function [outStats] = getClusterAnalysis(stats,p,d,collect, data, clusters, estimatedLabels, reducedParSpace)
%GETCLUSTERANALYSIS Summary of this function goes here
%   Detailed explanation goes here

[sortedFitness,id] = sort(collect.raw.pred.fitness_true(:,end,:));
if length(data)>1
    % idmap contains best run IDs for every bin
    idmap = id(1,:);
    % Do not compare when all models produced invalid shape in bin
    idmap(logical(squeeze(sum(collect.raw.pred.invalid(:,end,:)) == length(stats.cfg.xpID) ))) = NaN;
    clusterIDMap = nan(625,1);
    for i=1:max(idmap(:))
        tMap = nan(625,1);
        tMap(~stats.agg.nanSolutions{i}) = estimatedLabels(stats.agg.optimaIDs{i});
        clusterIDMap(~stats.agg.nanSolutions{i}.*idmap==i) = tMap(~stats.agg.nanSolutions{i}.*idmap==i);
    end
else
    tMap = nan(625,1);
    tMap(~stats.agg.nanSolutions{1}) = estimatedLabels(stats.agg.optimaIDs{1});
    clusterIDMap(~stats.agg.nanSolutions{1}) = tMap(~stats.agg.nanSolutions{1});
end


allNans = cell2mat(stats.agg.nanSolutions);



% For every binID
binIDs = find(ismember(clusterIDMap,clusters.uniqid));
largeClusterBinIDs = clusterIDMap;
largeClusterBinIDs(~ismember(clusterIDMap,clusters.uniqid)) = 0;
largeClusterBinIDs(ismember(clusterIDMap,clusters.uniqid)) = binIDs;
reducedParSpaceWithNaNs = nan(length(allNans),size(reducedParSpace,2)); reducedParSpaceWithNaNs = double(reducedParSpaceWithNaNs);
reducedParSpaceWithNaNs(allNans==0,:) = reducedParSpace;

bestOptimaLocations = [];
for binID=1:625
    %   Check segmentID
    segmentID = clusterIDMap(binID);
    if ismember(segmentID,clusters.uniqid)
        %   Get sample
        if length(data)>1
            %       Get runID
            runID = idmap(binID);
        else
            runID = 1;
        end
        %       Get sample location from run
        bestOptimaLocations(binID,:) = reducedParSpaceWithNaNs((runID-1)*625+binID,:);
    end
end


samplesLinIndices = stats.samples.mapLinIndx';samplesLinIndices=samplesLinIndices(:);
sampleFitnesses = [];
for m=1:length(stats.cfg.xpID); sampleFitnesses = [sampleFitnesses data{m}.samples.fitness_true]; end

% Generate output
outStats.sortedFitness = sortedFitness;
%outStats.idmap = idmap;
%outStats.clusterIDMap = clusterIDMap;
outStats.allNans = allNans;
%outStats.binIDs = binIDs;
outStats.largeClusterBinIDs = largeClusterBinIDs;
%outStats.reducedParSpaceWithNaNs =reducedParSpaceWithNaNs;
outStats.bestOptimaLocations = bestOptimaLocations;
outStats.samplesLinIndices = samplesLinIndices;
outStats.sampleFitnesses = sampleFitnesses;
end

