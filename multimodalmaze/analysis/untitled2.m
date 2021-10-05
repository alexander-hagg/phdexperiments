set(0,'DefaultFigureWindowStyle','default');
runID = 1;

trajectories = runData{runID}.acqMap.phenotype;
v1 = reshape(runData{runID}.acqMap.alignment,900,1);
v2 = reshape(runData{runID}.acqMap.traject2Goal,900,1);
v3 = reshape(runData{runID}.acqMap.firstGoal,900,1);
values = [v1(any(~isnan(v1),2)) v2(any(~isnan(v2),2)) v3(any(~isnan(v3),2))];
values = v3(any(~isnan(v3),2));
figHandles = vClassTrajectories(trajectories,values,runData{runID}.classification.labels,runData{runID}.d,runData{runID}.selection,3)


%%



save_figures(figHandles, '.', 'classes', 12, [5 4])

