%%
folder = '/scratch/ahagg2s/GECCO2019/planner_MUTATION_4096';
%folder = '/scratch/ahagg2s/GECCO2019/controller_MUTATION_4096';
analysisMutrate(folder);
load([folder '/analysis.mat']);

fig(1) = figure(1); hold off; fig(2) = figure(2); hold off; fig(3) = figure(3); hold off;
clrs = parula(4); style = {'--',':','-'};


%% Significance tests
clear signif;
sigLevels = [0.05, 0.005];
sigLevelID = 1;
%for sigLevelID = 1:length(sigLevels)
    sigLevel = sigLevels(sigLevelID);
    for mutID = 1:8
        signif(1,mutID) = kstest2(squeeze(exitCorrectness(3,mutID,:)),squeeze(exitCorrectness(2,mutID,:)),'Alpha',sigLevel);
        signif(2,mutID) = kstest2(squeeze(exitInCorrectness(3,mutID,:)),squeeze(exitInCorrectness(2,mutID,:)),'Alpha',sigLevel);
        signif(3,mutID) = kstest2(squeeze(userSelDrift(3,mutID,:)),squeeze(userSelDrift(2,mutID,:)),'Alpha',sigLevel);
    end
%end
signif


%% Mean and STD
median(exitCorrectness(:,4:end,:),[2:3])/900
median(exitInCorrectness(:,4:end,:),[2:3])/900
median(userSelDrift(:,4:end,:),[2:3])


%%
for methodID=1:length(methods)
    fig(1) = figure(1);
    pl(1,methodID) = semilogx(mutationRates,mCor(methodID,:),style{methodID},'Color',clrs(methodID,:),'LineWidth',3); hold on;
    scatter(mutationRates,mCor25(methodID,:),64,clrs(methodID,:),'filled','v'); hold on;
    scatter(mutationRates,mCor75(methodID,:),64,clrs(methodID,:),'filled','^'); hold on;
    ylabel('% correct');
    xlabel('Mutation Rate');xticks(mutationRates);xticklabels(mutationRates);grid on;axis([mutationRates(1) mutationRates(end)+0.01 0 100]);
    
    fig(2) = figure(2);
    pl(2,methodID) = semilogx(mutationRates,mInCor(methodID,:),style{methodID},'Color',clrs(methodID,:),'LineWidth',3); hold on;
    xlabel('Mutation Rate');xticks(mutationRates);xticklabels(mutationRates);grid on;axis([mutationRates(1) mutationRates(end)+0.01 0 100]);
    ylabel('% Incorrect');
    scatter(mutationRates,mInCor25(methodID,:),64,clrs(methodID,:),'filled','v'); hold on;
    scatter(mutationRates,mInCor75(methodID,:),64,clrs(methodID,:),'filled','^'); hold on;
    
    fig(3) = figure(3);
    pl(3,methodID) = loglog(mutationRates,mDrift(methodID,:),style{methodID},'Color',clrs(methodID,:),'LineWidth',3); hold on;
    xlabel('Mutation Rate');xticks(mutationRates);xticklabels(mutationRates);grid on;axis([mutationRates(1) mutationRates(end)+0.01 0 1]);
    ylabel('User Selection Drift');
    scatter(mutationRates,mDrift25(methodID,:),64,clrs(methodID,:),'filled','v'); hold on;
    scatter(mutationRates,mDrift75(methodID,:),64,clrs(methodID,:),'filled','^'); hold on;
    
    drawnow;
    
end


figure(1);legend(pl(1,:),methods,'Location','SouthEast');
figure(2);legend(pl(2,:),methods);
figure(3);legend(pl(3,:),methods,'Location','SouthEast');
drawnow;

%%
save_figures(fig, '.', 'mutRateComparison_Controller', 18, [8 5]);






