%% Show trajectories in allMaps
for selectionID=1%:size(allMaps,1)
    for runID=[1 size(allMaps,2)]
        maps = squeeze([allMaps(selectionID,runID,:)]);empty = cellfun(@isempty,maps);
        lastMap = find(empty); if ~isempty(lastMap); lastMap = lastMap(1)-1;else;lastMap=length(maps);end
        for mapID=lastMap
            genes = maps{lastMap}.genes; genes = reshape(genes,size(genes,1)*size(genes,2),size(genes,3)); genes = genes(~any(isnan(genes')),:);
            [trajectories,values] = eval_maze(genes,d{runID}.numHidden,d{runID}.maze,d{runID}.timesteps,d{runID}.ncores,d{runID}.goalLocations,d{runID}.useRNN,d{runID}.debug,false,true);    
            vTrajectories(trajectories, values{3}, d{runID}.maze, true)
        end
    end
end