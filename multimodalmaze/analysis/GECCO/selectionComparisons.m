set(0,'DefaultFigureWindowStyle','docked');
d = domain_Maze;
angles = [30 90 150 210 270 330];
numA =length(angles);

offset = 0;
offset = 3
figID = 1;
for i=1:length(data)/numA
    for angle=2%1:numA
        figure(figID+offset);
        dataO = data((angle-1)*length(data)/numA + i);
        subplot(1,3,1);
        viewMap(dataO.maps{1}.ring3,d);
        subplot(1,3,2);
        viewMap(dataO.maps{2}.ring3,d);
        subplot(1,3,3);
        viewMap(dataO.maps{end}.ring3,d);
        title(['Exit ' int2str(dataO.description.ringID) ', angle ' int2str(dataO.description.angle)]);
        colormap(parula(4));
        drawnow;
        figID = figID + 1;
        
    end
end

%%
set(0,'DefaultFigureWindowStyle','default');
d = domain_Maze;
angles = [30 90 150 210 270 330];
numA =length(angles);

clear summedExits exitCorrectness meanVals prctileVals prctileNonVals;
exits = [1 2 3];
experimentDirs = {'/scratch/ahagg2s/GECCO2019/planner_UDHM/selectedMaps.mat','/scratch/ahagg2s/GECCO2019/planner_SEEDS/selectedMaps.mat'};
for eD=1:length(experimentDirs)
    load(experimentDirs{eD});
    data = data(:,:,angles);
    % data = data(:,~cellfun(@isempty,{data(1,:).fname}));
    
    for exitID=[2:4]
        for angle=1:numA
            dataO = data(exitID,angle);
            if isempty(dataO.fname); continue; end
            for mapID=1:length(dataO.maps(end,:))
                exitCorrectness(eD,mapID,angle,1) = sum(dataO.maps{end,mapID}.ring1(:)==dataO.description.ringID);
                exitCorrectness(eD,mapID,angle,2) = sum(ismember( dataO.maps{end,mapID}.ring1(:), setdiff(exits,dataO.description.ringID)) );
            end
        end
    end
end

% subselection
subselectGens = [1 6 10];
exitCorrectness = exitCorrectness(:,subselectGens,:,:);


for s=1:size(exitCorrectness,1)
    meanVals(s,:,1) = mean(squeeze(exitCorrectness(s,:,:,1)),2);
    meanVals(s,:,2) = mean(squeeze(exitCorrectness(s,:,:,2)),2);
    prctileVals(s,:,1) = prctile(squeeze(exitCorrectness(s,:,:,1)),25,2);
    prctileVals(s,:,2) = prctile(squeeze(exitCorrectness(s,:,:,1)),75,2);
    prctileNonVals(s,:,1) = prctile(squeeze(exitCorrectness(s,:,:,2)),25,2);
    prctileNonVals(s,:,2) = prctile(squeeze(exitCorrectness(s,:,:,2)),75,2);
end

percCorrectVals = meanVals(:,:,1)./(meanVals(:,:,1) + meanVals(:,:,2));
percprctileCorrect25 = prctileVals(:,:,1)./(meanVals(:,:,1) + meanVals(:,:,2));
percprctileCorrect75 = prctileVals(:,:,2)./(meanVals(:,:,1) + meanVals(:,:,2));
%percprctileCorrect = [percprctileCorrect25;percprctileCorrect75];

numCorrectVals = meanVals(:,:,1);
numInCorrectVals = meanVals(:,:,2);
numprctileCorrect25 = prctileVals(:,:,1);
numprctileCorrect75 = prctileVals(:,:,2);
numprctileCorrect = [numprctileCorrect25;numprctileCorrect75];

%
lw = 2;style = {'k-','k--','k:'};
fig(1) = figure(1);hold off;
for s=1:size(exitCorrectness,1)
    pl(1,s) = plot(100*numCorrectVals(s,:,:)./900',style{s},'LineWidth',lw); hold on;
    plot([1:size(prctileVals,2);1:size(prctileVals,2)]+(s-1)*0.05  , 100*[prctileVals(s,:,1);prctileVals(s,:,2)]./900,style{s},'LineWidth',lw); hold on;
    plot([1:size(prctileVals,2);1:size(prctileVals,2)]+(s-1)*0.05  , 100*[prctileVals(s,:,1);prctileVals(s,:,2)]./900,style{s},'LineWidth',lw); hold on;
    
    pl(2,s) = plot(100*numInCorrectVals(s,:,:)./900',style{s},'LineWidth',lw,'Color',[1 0 0]); hold on;
    plot([1:size(prctileNonVals,2);1:size(prctileNonVals,2)]+(s+1)*0.05  , 100*[prctileNonVals(s,:,1);prctileNonVals(s,:,2)]./900,style{s},'LineWidth',lw,'Color',[1 0 0]); hold on;
    plot([1:size(prctileNonVals,2);1:size(prctileNonVals,2)]+(s+1)*0.05  , 100*[prctileNonVals(s,:,1);prctileNonVals(s,:,2)]./900,style{s},'LineWidth',lw,'Color',[1 0 0]); hold on;
end
grid on;
ylabel('% Correct Exits');
xlabel('Generations');
axis([0.8 length(subselectGens)+0.2 0 100]);
ax = gca;
ax.XTick = 1:length(subselectGens);
ax.XTickLabel = floor(0:(8192/(length(subselectGens)-1)):8192);
legend([pl(1,:),pl(2,:)], 'UDHM S','Seed S','UDHM not S','Seed not S','Location','NorthWest');
drawnow;
save_figures(fig, '.', 'seed+udhm+combined', 24, [8 8]);


%%
angles = [30 90 150 210 270 330];

for exitID=1:3
    for angleID=1:length(angles)
        angle = angles(angleID);
        caseStr = int2str(exitID);
        dataPath = ['GECCO2019/runsSEED/exit' int2str(exitID) '/' int2str(angle) '.mat']
        load(dataPath);
        figure(1);
        viewMap(runData{end}.acqMap.ring3,runData{end}.d)
        title(['SEED: ' caseStr ' - selThresh ' num2str(runData{2}.selThresh)]);
        caxis([0 3]); colormap(parula(4));
        figure(2);
        viewMap(runData{end}.acqMap.fitness,runData{end}.d)
        colormap(parula(10));
        
        dataPath = ['GECCO2019/runsUDHM/exit' int2str(exitID) '/' int2str(angle) '.mat']
        load(dataPath);
        figure(3);
        viewMap(runData{end}.acqMap.ring3,runData{end}.d)
        caxis([0 3]); colormap(parula(4));
        title(['UDHM: ' caseStr ' - selThresh ' num2str(runData{2}.selThresh)]);
        figure(4);
        viewMap(runData{end}.acqMap.fitness,runData{end}.d)
        colormap(parula(10));
        drawnow;
        pause(2);
    end
end


%%
samples = reshape(runData{end}.acqMap.genes,900,92);
[adjustedFitness, trajectories, values] = runData{end}.d.fitfun(samples(~any(isnan(samples')),:));
%%
figure(10)
vTrajectories(trajectories, values{4}, runData{1}.d.mazeFileName, true);


