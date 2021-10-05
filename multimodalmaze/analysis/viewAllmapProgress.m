set(0,'DefaultFigureWindowStyle','docked');
maxFitVal = 500;
showCuriousness = false;

for selectionID=1%:size(allMaps,1)
    for runID=[size(allMaps,2)-1 size(allMaps,2)]
        maps = squeeze([allMaps(selectionID,runID,:)]);empty = cellfun(@isempty,maps);
        lastMap = find(empty); if ~isempty(lastMap); lastMap = lastMap(1)-1;else;lastMap=length(maps);end
        for mapID=1:lastMap 
            if (lastMap ==0) || isempty(allMaps{selectionID,runID,mapID});break;end
            fig((runID-1)*3 + 1) = figure((runID-1)*3 + 1);
            fitness = allMaps{selectionID,runID,mapID}.fitness;
            [~,~,cH] = viewMap(fitness,d{runID});
            cH.Label.String = 'Fitness';
            caxis([0 maxFitVal]);
            title(['Run ID: ' int2str(runID) ' - map ' int2str(mapID)]);
            
            if showCuriousness
                fig((runID-1)*3 + 2) = figure((runID-1)*3 + 2);
                curiousness = allMaps{selectionID,runID,mapID}.curiousness;
                curiousness(isnan(fitness)) = nan;
                [~,~,cH] = viewMap(curiousness,d{runID});
                cH.Label.String = 'Curiousness';
                %caxis([-30 0]);
            end
            
            fig((runID-1)*3 + 3) = figure((runID-1)*3 + 3);
            firstGoal = allMaps{selectionID,runID,mapID}.firstGoal;
            firstGoal(isnan(fitness)) = nan;
            [~,~,cH] = viewMap(firstGoal,d{runID});
            cH.Label.String = 'Exit ID';
            caxis([-0.5 3.5]);colormap(parula(4));
            cH.Ticks = [0 1 2 3];cH.TickLabels = {'none','1','2','3'};
            if selectionID < 4
                title(['Selected Exit: ' int2str(selectionID)]);
            else
                title(['Selected Exit: none']);
            end
            %drawnow;
            %pause(0.1)
        end
    end
end
%%
save_figures(fig, '.', ['selectionOnMAP_ID_' int2str(selectionID)], 12, [5 4])
