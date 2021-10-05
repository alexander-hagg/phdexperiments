figure(1);hold off;
clrs = [1 0 0; 0 1 0; 0 0 1];
lw = [1 4 8];
p.minSelectedSamples = 50;
p.selectionThreshold = 1.01;
thresholds = 1:-0.01:0.01;
for runID = 2:3
    selStat = [];
    classification = runData{runID}.classification;
    p = runData{runID}.p;
    d = runData{runID}.d;
    values{1} = runData{runID}.acqMap.alignment;
    values{2} = runData{runID}.acqMap.traject2Goal;
    values{3} = runData{runID}.acqMap.firstGoal;
    for selVal =1:3
        for selThresID = 1:length(thresholds)
            selThresh = thresholds(selThresID);
            [isSelected,selection] = mazeSelection(classification.labels, values, ...
                p.selectionCriterionID, selVal, selThresh);
            selStat(selThresID,1) = selThresh;
            selStat(selThresID,2) = sum(isSelected);
            if (selThresh < p.selectionThreshold) && sum(isSelected) > p.minSelectedSamples
                disp(['Threshold at: ' num2str(selThresh) ', selected: ' int2str(sum(isSelected))]);
                %break;
            end
        end
        
        h(selVal) = plot(selStat(:,2),selStat(:,1),'Color',clrs(selVal,:),'LineWidth',lw(runID));
        axis([0 800 0 1]);
        grid on; grid minor;
        hold on;
    end
    legend([h(1) h(2) h(3)],'exit 1', 'exit 2', 'exit 3');
end
title(p.selectionValue);
