set(0,'DefaultFigureWindowStyle','docked');


exitStrs = {'1','2','3','none'};
xtickpos = [0 2^12];
xtickLabels = {'2^{16}','2^{16}+2^{12}'};
maxVal = 500;

selectionGens = [1];

firstGoalSeq =[];
for j=1:size(firstGoals,2)
    for i=1:size(firstGoals,1)
        firstGoalSeq(1,i,j,:) = firstGoals{i,j}(:,1);
        firstGoalSeq(2,i,j,:) = firstGoals{i,j}(:,2);
        firstGoalSeq(3,i,j,:) = firstGoals{i,j}(:,3);
    end
end


%
clrs = parula(4);
clrs(1,:) = [];
for j=1:size(firstGoals,2)
    fig(j) = figure(j);hold off;
    mV = squeeze(firstGoalSeq(1,:,j,:));
    x2 = [[1:length(mean(mV))], fliplr([1:length(mean(mV))])];inBetween = [mean(mV)-std(mV), fliplr(mean(mV)+std(mV))];
    h = fill(x2', inBetween', clrs(1,:));hold on;set(h,'facealpha',.3);
    plI(1) = plot(mean(mV),'Color',clrs(1,:),'LineWidth',4);
    ax = gca;
    ax.XTick = xtickpos;
    ax.XTickLabel = xtickLabels;
    
    %
    plot(mean(squeeze(firstGoalSeq(2,:,j,:))),'Color',clrs(2,:),'LineWidth',2);hold on;
    mV = squeeze(firstGoalSeq(2,:,j,:));
    x2 = [[1:length(mean(mV))], fliplr([1:length(mean(mV))])];inBetween = [mean(mV)-std(mV), fliplr(mean(mV)+std(mV))];
    h = fill(x2', inBetween', clrs(2,:));hold on;set(h,'facealpha',.3);
    plI(2) = plot(mean(mV),'Color',clrs(2,:),'LineWidth',4);
    
    plot(mean(squeeze(firstGoalSeq(3,:,j,:))),'Color',clrs(3,:),'LineWidth',2);hold on;
    mV = squeeze(firstGoalSeq(3,:,j,:));
    x2 = [[1:length(mean(mV))], fliplr([1:length(mean(mV))])];inBetween = [mean(mV)-std(mV), fliplr(mean(mV)+std(mV))];
    h = fill(x2', inBetween', clrs(3,:));hold on;set(h,'facealpha',.2);
    plI(3) = plot(mean(mV),'Color',clrs(3,:),'LineWidth',4);

    plII = plot([selectionGens' selectionGens']',[0 maxVal]','k--');
    grid on;axis([-10 size(firstGoalSeq,4) 0 maxVal]);
    grid minor;
    title(['Exit selected ' exitStrs{j}]);
    xlabel('Generations'); ylabel('# Elites');
    legend([plI(1) plI(2) plI(3) plII(1)],'Exit 1','Exit 2','Exit 3','Selection','Location','NorthWest');
end


%%
save_figures(fig, '.', 'selectionExample_Hid5_Sobol_Sobol', 18, [5 4])
