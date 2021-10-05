clear;clc;

%domainname = 'PARSEC';systemInit;
%[names,data] = read_experiments('/scratch/ahagg2s/acq_parsec2/1');

domainname = 'FOILFFD';systemInit;
[names,data] = read_experiments('/scratch/ahagg2s/acq_gps_long5+6/1');


%% Analysis
d = data{1}.d;p = data{1}.p;
disp(['Retrieving statistics']);
[stats,p,d,collect] = getModelComparisonStats(data, p, d, names, 'bindiff', 0.03);


%% Visualize comparison between model types
set(0,'DefaultFigureWindowStyle','docked')
colors = {'r','b', 'g', 'k'};
uniqueXPs = 1:length(stats.cfg.xpUNIQ);

for lo = uniqueXPs(3)
    leaveout = lo;
    selmaps = setxor(uniqueXPs,leaveout);
    for j = 1:length(selmaps)
        i = selmaps(j);
        [sortedRelativeFitness, sortedFitnessIDs] = sort(squeeze(stats.pred.bins.relativeFitness(i,lo,:)));
        
        fig((j-1)*4 + 1) = figure((j-1)*4 + 1);hold off;
        diffmap = nan(d.featureRes(1),d.featureRes(2));
        diffmap(~stats.pred.bins.invalid(i,lo,:)) = 0;
        diffmap(:) = diffmap(:) + squeeze(stats.pred.bins.better(i,lo,:));
        diffmap(:) = diffmap(:) - squeeze(stats.pred.bins.worse(i,lo,:));
        [f,~,h] = viewMap(diffmap, d);
        caxis([-1 1]);
        h.Ticks = [-0.7 0 0.7];
        h.TickLabels = {['< ' num2str(100*-stats.cfg.bindiff) '% (' int2str(stats.summary.bins.worse(i,lo)) ')'], 'similar' , ['> ' num2str(100*stats.cfg.bindiff) '% (' int2str(stats.summary.bins.better(i,lo)) ')']};
        colormap([[1.0 0 0]; [0.5 0.5 0.5]; [0.0 1.0 0.0]]);
        t = title(['Fitness Comparison ' stats.cfg.xpUNIQ{i} ' vs ' stats.cfg.xpUNIQ{lo}]);
        
        fig((j-1)*4 + 2) = figure((j-1)*4 + 2);hold off;
        densitydifferencemap = squeeze(stats.samples.meanMapLocations(i,:,:)) - squeeze(stats.samples.meanMapLocations(lo,:,:));
        [f,~,h] = viewMap(densitydifferencemap, d);
        caxis([-2 2]);
        h.Ticks = [-1 0 1];
        h.TickLabels = {stats.cfg.xpUNIQ{lo},'equal',stats.cfg.xpUNIQ{i}};
        %colormap([[1.0 0 0]; [1 0.5 0.5] ; [0.5 0.5 0.5]; [0.5 1 0.5] ; [0.0 1.0 0.0]]);
        colormap(flipud(redgreencmap(50, 'Interpolation', 'quadratic')));
        title(['Difference in Sample Density']);
        
        
        fig((j-1)*4 + 3) = figure((j-1)*4 + 3);hold off;
        [f,~,h] = viewMap(diffmap.*(densitydifferencemap), d);
        caxis([-1 1]);
        %h.Ticks = [-0.7 0 0.7];
        %h.TickLabels = {['< ' num2str(100*-stats.cfg.bindiff) '% (' int2str(stats.summary.bins.worse(i,lo)) ')'], 'similar' , ['> ' num2str(100*stats.cfg.bindiff) '% (' int2str(stats.summary.bins.better(i,lo)) ')']};
        %colormap([[1.0 0 0]; [0.5 0.5 0.5]; [0.0 1.0 0.0]]);
        t = title(['Density*Fitness ' stats.cfg.xpUNIQ{i} ' vs ' stats.cfg.xpUNIQ{lo}]);
        
        
        fig((j-1)*4 + 4) = figure((j-1)*4 + 4);hold off;
        colorscatter = [repmat([1.0 0.0 0.0],sum(sortedRelativeFitness<-stats.cfg.bindiff),1) ; ...
            repmat([0.0 0.0 0.0],sum(sortedRelativeFitness<stats.cfg.bindiff)-sum(sortedRelativeFitness<-stats.cfg.bindiff),1) ; ...
            repmat([0.0 1.0 0.0],625-sum(sortedRelativeFitness<stats.cfg.bindiff),1)];
        scatter(1:length(sortedRelativeFitness),sortedRelativeFitness, [], colorscatter);
        hold on;
        plot([0 625], [0 0],'k:', 'LineWidth', 2);
        plot([0 625], [-stats.cfg.bindiff -stats.cfg.bindiff],'k:', 'LineWidth', 2);
        plot([0 625], [stats.cfg.bindiff stats.cfg.bindiff],'k:', 'LineWidth', 2);
        axis([0 625 -3 3]);
        ax = gca;
        ax.YTick = [-stats.cfg.bindiff 0 stats.cfg.bindiff];
        ax.YTickLabel = string(100*ax.YTick);
        grid on;
        xlabel('Bins');
        ylabel(['% higher bin fitness']);
        title([stats.cfg.xpUNIQ{i} ' vs ' stats.cfg.xpUNIQ{lo}]);
        hold off;
        
        
        drawnow;
    end
    %pause(1)
    %save_figures(fig, './', ['binfitAndParameters_' int2str(lo) '_'], 12, [5 5]);
    
end

%% Visualize comparison between model instances
figure(1);hold off;
for i=1:8
    viewMap(squeeze(stats.samples.mapLocations(i,:,:)),d);
    caxis([0 10]);
    pause(1);
end
%% Compare ranges and max, min
clrs = parula(length(stats.cfg.xpUNIQ));
fade = 1:-1/(length(stats.cfg.xpID)/length(stats.cfg.xpUNIQ)):1/(length(stats.cfg.xpID)/length(stats.cfg.xpUNIQ));
fade = repmat(fade,1,length(stats.cfg.xpUNIQ));
clrs = repelem(clrs,length(stats.cfg.xpUNIQ),1);
clrs = clrs.*fade';
names = [stats.cfg.xpUNIQ, 'All'];
for i=1:length(names)
    if i < length(stats.cfg.xpUNIQ)+1
        groups = [stats.cfg.xpID==i];
    else
        groups = [stats.cfg.xpID~=0];
    end
    figure((i-1)*4 + 1)
    [f,~,h] = viewMap(100*reshape(squeeze(nanmax(collect.raw.pred.fitness_true_rel(groups,end,:))),25,25)./ ...
        reshape(squeeze(nanmax(collect.raw.pred.fitness_true_rel([stats.cfg.xpID~=0],end,:))),25,25),d);
    h.Label.String = '% of best found in bin';
    caxis([90 100]);
    colormap(flipud(redgreencmap(10, 'Interpolation', 'linear')));
    title(['Best rel. Fitness (' names{i} ')']);
    
    figure((i-1)*4 + 2)
    [f,~,h] = viewMap(100*reshape(squeeze(nanmin(collect.raw.pred.fitness_true_rel(groups,end,:))),25,25)./ ...
        reshape(squeeze(nanmax(collect.raw.pred.fitness_true_rel([stats.cfg.xpID~=0],end,:))),25,25),d);
    h.Label.String = '% of best found in bin';
    caxis([80 100]);
    colormap(flipud(redgreencmap(10, 'Interpolation', 'linear')));
    title(['Worst rel. Fitness (' names{i} ')']);
    
    figure((i-1)*4 + 3)
    [f,~,h] = viewMap(reshape(squeeze(range(collect.raw.pred.fitness_true_rel(groups,end,:))),25,25),d);
    h.Label.String = '% of best found in bin';
    title(['Fitness % Range (' names{i} ')']);
    caxis([0 0.3]);
    colormap(redgreencmap(10, 'Interpolation', 'linear'));
    
    [~,id] = sort(collect.raw.pred.fitness_true(groups,end,:));
    id = squeeze(id);
    figure((i-1)*4 + 4)
    idmap = id(1,:);
    % Do not compare when all models produced invalid shape in bin
    idmap(logical(squeeze(sum(collect.raw.pred.invalid(groups,end,:)) == sum(groups) ))) = NaN;
    [f,~,h] = viewMap(reshape(idmap,25,25),d);
    title(['Which model was best?']);
    caxis([1 sum(groups)]);
    h.Ticks = 1:sum(groups);
    if i < length(names)
        labels{i} = strcat([names{i} '_'], h.TickLabels);
        h.TickLabels = labels{i};
    else
        h.TickLabels = cat(1,labels{:});
    end
    colormap(clrs(groups,:));
    drawnow;
end


%%
fig(1) = figure(101);hold off;
%
cmap = colorcube(ceil(max(estimatedLabels))+1);clrs = cmap(estimatedLabels+1,:);clrs(estimatedLabels==0,:) = repmat([0 0 0],sum(estimatedLabels==0),1);
for i=1:length(selected)
    disp(['Method found ' num2str(100*sum(unique(estimatedLabels(pickedOptima{i})))/sum(unique(estimatedLabels))) '% of total clusters']);
    scatter(reducedParSpace(pickedOptima{i},1),reducedParSpace(pickedOptima{i},2),[],clrs(pickedOptima{i},:),'filled');
    hold on;axis equal;axis tight;
end

ax = gca;ax.XTickLabel = [];ax.YTickLabel = [];ax.ZTickLabel = [];
title(['Optima (t-SNE dimensionality reduction, DBSCAN clustering)']);

save_figures(fig, './', ['distribution_optima'], 14, [10 10]);