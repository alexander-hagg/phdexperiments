% Build random Gaussian dataset
clear;clc;
numSets = 100;
perplexities = [5 10 20 50 100];
dimensionalities = [2 5 10 20];

for d=1:length(dimensionalities)
    dims = dimensionalities(d)
    setPts = 3+randi(17,numSets,1);
    input = []; rlLab = [];
    for pp=1:length(setPts)
        width = 3*rand(1);
        input = [input; randi(10,1,dims)-5 + width*randn(setPts(pp),dims)];
        rlLab = [rlLab; pp*ones(setPts(pp),1)];
    end
    in{d} = input;
    realLabels{d} = rlLab;
    
    for i=1:length(perplexities)
        perplexity = perplexities(i)
        
        disp('Reduce Dimensionality');
        cfg.no_dims = 2; cfg.initial_dims = dims; cfg.theta = 0.5; cfg.alg = 'svd'; cfg.max_iter = 1000;
        reducedParSpace{d,i} = fast_tsne(in{d}, cfg.no_dims, cfg.initial_dims, perplexity, cfg.theta, cfg.alg, cfg.max_iter);
        
        coreneighbours = max(2*dims,3); %Rule of thumb
        [~,t_distances] = knnsearch(in{d},in{d},'K',coreneighbours+1);
        t_distances(:,1) = [];
        distances{1,d,i} = sort(t_distances(:));
        [maxVal{1,d,i},maxID{1,d,i}] = getElbow(distances{1,d,i});
        epsilon(1,d,i) = maxVal{1,d,i};
        tic;    [~,estimatedLabels{1,d,i},cen] = dbscan(in{d}', epsilon(1,d,i), coreneighbours); timing{1,d,i} = toc;
        randIndex(1,d,i) = rand_index(realLabels{d}, estimatedLabels{1,d,i}, 'adjusted');
        
        coreneighbours = max(2*2,3); %Rule of thumb
        [~,t_distances] = knnsearch(reducedParSpace{d,i},reducedParSpace{d,i},'K',coreneighbours+1);
        t_distances(:,1) = [];
        distances{2,d,i} = sort(t_distances(:));
        [maxVal{2,d,i},maxID{2,d,i}] = getElbow(distances{2,d,i});
        epsilon(2,d,i) = maxVal{2,d,i};
        tic;    [~,estimatedLabels{2,d,i},cen] = dbscan(reducedParSpace{d,i}', epsilon(2,d,i), coreneighbours); timing{2,d,i} = toc;
        randIndex(2,d,i) = rand_index(realLabels{d}, estimatedLabels{2,d,i}, 'adjusted');
        
    end
end

epsilon
randIndex
%% Visualization
set(0,'DefaultFigureWindowStyle','default');

fig(1) = figure(1);hold off;
subplot(2,1,1);
plot(squeeze(randIndex(1,:,:)), 'LineWidth',2)

legend(strcat('Perp: ', string(perplexities)), 'Location', 'NorthWest');
ax = gca;ax.XTick = 1:max(ax.XTick);
ax.XTickLabel = string(dimensionalities);grid on;xlabel('Dimensionality');ylabel('Rand Index');
title('DBSCAN>REDUX');

subplot(2,1,2);
plot(squeeze(randIndex(2,:,:)), 'LineWidth',2)

legend(strcat('Perp: ', string(perplexities)), 'Location', 'NorthWest');
ax = gca;ax.XTick = 1:max(ax.XTick);
ax.XTickLabel = string(dimensionalities);grid on;xlabel('Dimensionality');ylabel('Rand Index');
title('REDUX>DBSCAN');



for d=1:length(dimensionalities)
    
    for i=1:length(perplexities)
        fig(2) = figure(2);hold off;
        plot(distances{1,d,i});hold on;
        scatter(maxID{1,d,i},epsilon(1,d,i),'filled');
        colormap(hsv);
        legend('original','max','Location','NorthWest');ylabel('Pair-Wise Distance');xlabel('Ascending Order');
        grid on;title('NN distances original space');
        
        fig(3) = figure(3);hold off;
        plot(distances{2,d,i});hold on;
        scatter(maxID{2,d,i},epsilon(2,d,i),'filled');
        colormap(hsv);
        legend('original','max','Location','NorthWest');ylabel('Pair-Wise Distance');xlabel('Ascending Order');
        grid on;title('NN distances reduced space');
        
        
        clrsREAL = hsv(max(realLabels{d}));
        clrs = hsv(max(estimatedLabels{1,d,i})+1);
        clrs(1,:) = zeros(3,1);
        clrsREDUCED = hsv(max(estimatedLabels{2,d,i})+1);
        clrsREDUCED(1,:) = zeros(3,1);
        dotsize = 8;
        
        fig(4) = figure(4);hold off;
        h1 = scatter(squeeze(in{d}(:,1))', squeeze(in{d}(:,2))',dotsize,clrsREAL(realLabels{d},:),'filled');
        title([int2str(numSets) ' Gaussians clusters']);
        grid on;axis equal;
       
        fig(5) = figure(5);hold off;
        h1 = scatter(squeeze(reducedParSpace{d,i}(:,1))', squeeze(reducedParSpace{d,i}(:,2))',dotsize,clrsREAL(realLabels{d},:),'filled');
        title([int2str(numSets) ' Gaussians clusters']);
        grid on;axis equal;
        
        fig(6) = figure(6);hold off;
        h1 = scatter(squeeze(reducedParSpace{d,i}(estimatedLabels{1,d,i}==0,1))', squeeze(reducedParSpace{d,i}(estimatedLabels{1,d,i}==0,2))',dotsize,clrs(estimatedLabels{1,d,i}(estimatedLabels{1,d,i}==0)+1,:),'filled');
        hold on;
        h2 = scatter(squeeze(reducedParSpace{d,i}(estimatedLabels{1,d,i}>0,1))', squeeze(reducedParSpace{d,i}(estimatedLabels{1,d,i}>0,2))',dotsize,clrs(estimatedLabels{1,d,i}(estimatedLabels{1,d,i}>0)+1,:),'filled');
        title(['DBSCAN>REDUX: ' int2str(max(estimatedLabels{1,d,i})) ' clusters, Rand index: ' num2str(randIndex(1,d,i))]);
        legend([h1 h2], [num2str(round(100*sum(estimatedLabels{1,d,i}==0)/length(estimatedLabels{1,d,i}))) '% unassigned'], [int2str(max(estimatedLabels{1,d,i})) ' clusters']);
        grid on;axis equal;
        
        fig(7) = figure(7);hold off;
        h1 = scatter(squeeze(reducedParSpace{d,i}(estimatedLabels{2,d,i}==0,1))', squeeze(reducedParSpace{d,i}(estimatedLabels{2,d,i}==0,2))',dotsize,clrsREDUCED(estimatedLabels{2,d,i}(estimatedLabels{2,d,i}==0)+1,:),'filled');
        hold on;
        h2 = scatter(squeeze(reducedParSpace{d,i}(estimatedLabels{2,d,i}>0,1))', squeeze(reducedParSpace{d,i}(estimatedLabels{2,d,i}>0,2))',dotsize,clrsREDUCED(estimatedLabels{2,d,i}(estimatedLabels{2,d,i}>0)+1,:),'filled');
        title(['REDUX>DBSCAN: ' int2str(max(estimatedLabels{2,d,i})) ' clusters, Rand index: ' num2str(randIndex(2,d,i))]);
        legend([h1 h2], [num2str(round(100*sum(estimatedLabels{2,d,i}==0)/length(estimatedLabels{2,d,i}))) '% unassigned'], [int2str(max(estimatedLabels{2,d,i})) ' clusters']);
        grid on;axis equal;
        drawnow;
    end
end

save_figures(fig, './', ['tsne_results_'], 12, [7 7]);