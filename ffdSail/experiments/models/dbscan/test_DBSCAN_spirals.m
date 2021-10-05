%data = twospirals(N, degrees, start, noise);
data = twospirals(2000);
input = data(:,1:2);
labels = data(:,3);

fig(1) = figure(1);
scatter(input(:,1),input(:,2),[],labels);
axis([-12 12 -12 12]);axis equal;grid on;
title('True Labels');

coreneighbours = max(2 * 2,3); 
[~,t_distances]     = knnsearch(input,input,'K',coreneighbours+1);
t_distances(:,1)    = [];
t_distances         = sort(t_distances(:));
[maxVal ,maxID]     = getElbow(t_distances);
%epsilon             = maxVal;

epsilonVals = [maxVal, 1, 2, 3];
clear centers estimations
for i=1:length(epsilonVals)
    epsilon = epsilonVals(i);
    tic; [~,estimatedLabels,cen] = dbscan(input', epsilon, coreneighbours); t = toc;
    estimations(i,:) = estimatedLabels';
end

for i=1:length(epsilonVals)
    fig(end+1) = figure(i+1);
    scatter(input(:,1),input(:,2),[],estimations(i,:));
    axis([-12 12 -12 12]);axis equal;grid on;
    title(['Estimated Labels, epsilon ' num2str(epsilonVals(i))]);
end

save_figures(fig, './', ['DBSCAN_spirals_'], 12, [5 5]);

