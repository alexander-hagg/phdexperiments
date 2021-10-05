set(0,'DefaultFigureWindowStyle','default');

for runID=1:1
    fig((runID-1)*2+1) = figure((runID-1)*2+1);
    [pltHandle] = viewClasses(runData{runID}.classification,true);
    axis([-10 10 -8 8]);
    axis equal;
    
    fig((runID-1)*2+2) = figure((runID-1)*2+2);
    
    [pltHandle] = viewClasses(runData{runID}.classification,false,runData{runID}.selection);
    axis([-10 10 -8 8]);
    axis equal;
    title('Selected');
    %%
    trajectories = runData{runID}.acqMap.phenotype;
    v1 = reshape(runData{runID}.acqMap.alignment,900,1);
    v2 = reshape(runData{runID}.acqMap.traject2Goal,900,1);
    v3 = reshape(runData{runID}.acqMap.firstGoal,900,1);
    values = [v1(any(~isnan(v1),2)) v2(any(~isnan(v2),2)) v3(any(~isnan(v3),2))];
    values = v3(any(~isnan(v3),2));
    figHandles = viewClassTrajectories(trajectories,values,runData{runID}.classification.labels,runData{runID}.d,runData{runID}.selection,3)
    
end



%%
save_figures(fig, '.', ['simspace_ID_' int2str(runData{runID}.p.selectionValue)], 12, [8 8])
%%
save_figures(figHandles, '.', 'classes', 12, [5 4])


