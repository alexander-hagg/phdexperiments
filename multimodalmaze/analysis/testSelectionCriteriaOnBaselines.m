% TODO Get Acquisition map from file
% folder = '/scratch/ahagg2s/PROPHESAI/ECJ_THETA_baselines/HID2_2k';
d = runData{1}.d; p = runData{1}.p;
map = runData{1}.acqMap;
classification = runData{1}.classification;
genes = reshape(map.genes,900,size(map.genes,3));genes(isnan(map.fitness(:)),:) = [];
%disp('Classification started'); 
%classification = extractClasses(genes); disp('Classification and prototyping done');


p.selectionThreshold = 0.8;

set(0,'DefaultFigureWindowStyle','docked');
values = [];
critID = 3;
firstGoal = [];
for i=1:size(map.phenotype,1)
    firstGoal(i) = evalFirstGoal(squeeze(map.phenotype(i,:,:)),d.goalLocations);
end
map.firstGoal(~isnan(map.firstGoal)) = firstGoal;

%
values{1} = map.alignment(~isnan(map.fitness(:)));
values{2} = map.traject2Goal(~isnan(map.fitness(:)));
values{3} = map.firstGoal(~isnan(map.fitness(:)));


minSelectedSamples = 10;
for selVal = 1:3
    disp('==================================')
    disp(selVal)
    thresholds = 1:-0.01:0.5;
    for selThresID = 1:length(thresholds)
        selThresh = thresholds(selThresID); 
        [isSelected,selection] = mazeSelection(classification.labels, values, ...
            critID, selVal, selThresh);        
        disp(['Selected ' int2str(sum(isSelected)) ' individuals.']);
        if (selThresh < p.selectionThreshold) && sum(isSelected) > minSelectedSamples
            disp(['Threshold at: ' num2str(selThresh)]);
            break;
        end
    end
    
    % Show selected invidiuals
    sel = ismember(classification.labels,selection);
    map.selected = nan(size(map.fitness));
    map.selected(~isnan(map.fitness)) = sel;
    fig(selVal) = figure(selVal);hold off;
    [~,~,cHandle] = viewMap(map.firstGoal.*map.selected,d);caxis([0 3]);
    title(['Selection Exit ' int2str(selVal)]);
    colormap(parula(4));
    cHandle.Ticks = 0:3;
    
    %fig(selVal+3) = figure(selVal+3);hold off;
    %constraints = setConstraints(classification, selection, p.constraintType);

    
    
    selectedGenesMap = reshape(map.genes,900,size(map.genes,3));
    selectedGenesMap = selectedGenesMap(~any(isnan(selectedGenesMap')),:);
    selectedGenesMap(~sel,:) = nan;
    sGP = reshape(map.genes,900,size(map.genes,3));
    sGP(~isnan(map.fitness),:) = selectedGenesMap;
    sGP(isnan(map.fitness),:) = nan;
    [trajectories] = d.evalFcn(sGP(~any(isnan(sGP')),:));
    
    fig(selVal+6) = figure(selVal+3);hold off;
    plotSim(trajectories,d,true,[0 0 0]);
    title(['Selection Exit ' int2str(selVal)]);
    drawnow;
end

%save_figures(fig, '.', 'selectionExamples', 12, [5 4])


