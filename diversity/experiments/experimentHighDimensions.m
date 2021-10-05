clear;clc;

dims = [1e0, 1e1, 1e2, 1e3, 1e4];
nSamples = 1000;

for i=1:numel(dims)
    dim = dims(i);
    samples = rand(nSamples,dim);
    distMetric = 'euclidean';
    distances = pdist2(samples,samples,distMetric);    
    fig(1) = figure(1);
    subplot(numel(dims),1,i);hold off;
    distances = distances./max(distances(:));
    histogram(distances(:),100);
    hold on;
    scatter(min(distances(:)),0,32,[1 0 0],'filled');
    scatter(max(distances(:)),0,32,[0 0 1],'filled');
    title(['Distances of ' int2str(nSamples) ' Random Samples, # Dims: ' int2str(dim)]);
    xlabel(['Distance: ' distMetric ' (normalized)']);
    ylabel('Count');
    legend('Distribution','Minimum','Maximum','Location','NorthWest');
    axis([0, 1, 0, 3e4 + 1e5*(i-1)]);
    drawnow;

    distMetric = 'cityblock';
    distances = pdist2(samples,samples,distMetric);    
    distances = distances./max(distances(:));
    fig(2) = figure(2);
    subplot(numel(dims),1,i);hold off;
    histogram(distances(:),100);
    hold on;
    scatter(min(distances(:)),0,32,[1 0 0],'filled');
    scatter(max(distances(:)),0,32,[0 0 1],'filled');
    title(['Distances of ' int2str(nSamples) ' Random Samples, # Dims: ' int2str(dim)]);
    xlabel(['Distance: ' distMetric ' (normalized)']);
    ylabel('Count');
    legend('Distribution','Minimum','Maximum','Location','NorthWest');
    axis([0, 1, 0, 3e4 + 1e5*(i-1)]);
    drawnow;

end

%%
save_figures(fig, '.', 'dimensionality', 10, [6 12])
