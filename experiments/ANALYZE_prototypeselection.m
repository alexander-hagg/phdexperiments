c = 1;
optima = reshape(data{1}.predMap(end/2).genes,625,10);
optima(any(isnan(optima')),:) = [];
    
clustermembers = optima(o{run}.estimatedLabels==clusters(c),:);
clusternonmembers = optima(o{run}.estimatedLabels~=clusters(c),:);

memberLocations = o{run}.reducedParSpace(o{run}.estimatedLabels==clusters(c),:);
nonmemberLocations = o{run}.reducedParSpace(o{run}.estimatedLabels~=clusters(c),:);


[a, ce] = kmedoids(memberLocations,1);


showdims = [5,6];

fig(1) = figure(1);hold off;
distances1 = pdist2(ce,memberLocations);
plot(distances1,'x');
hold on;
distances2 = pdist2(ce,nonmemberLocations);
plot(distances2,'o');
legend('members', 'non-members');
title('distances');

%
fig(2) = figure(2);hold off;
scatter(o{run}.reducedParSpace(:,1),o{run}.reducedParSpace(:,2),'filled');
hold on;
scatter(memberLocations(:,1),memberLocations(:,2),'filled');
scatter(ce(:,1),ce(:,2),128,'filled');
title('Cluster in reduced space');

%
fig(3) = figure(3);hold off;
scatter(optima(:,showdims(1)),optima(:,showdims(2)),'filled');
hold on;
scatter(clustermembers(:,showdims(1)),clustermembers(:,showdims(2)),'filled');
predmap = reshape(o{run}.sail{c}.predMap(end).genes,625,10);
scatter(predmap(:,showdims(1)),predmap(:,showdims(2)),'filled');
title('Cluster in 1st 2 dims of original space');
legend('500PE Optima', '500PE sel. cluster','1000PE PREDMAP');

%

fig(4) = figure(4);hold off;
scatter(optima(:,showdims(1)),optima(:,showdims(2)),'filled');
hold on;
scatter(clustermembers(:,showdims(1)),clustermembers(:,showdims(2)),'filled');
acqmap = reshape(o{run}.sail{c}.acqMap(end).genes,625,10);
scatter(acqmap(:,showdims(1)),acqmap(:,showdims(2)),'filled');
title('Cluster in 1st 2 dims of original space');
legend('500PE Optima', '500PE sel. cluster','1000PE ACQMAP');

%
% fig(5) = figure(5);hold off;
% scatter(o{run}.reducedParSpace(:,1),o{run}.reducedParSpace(:,2),'filled');
% hold on;
% scatter(memberLocations(:,1),memberLocations(:,2),'filled');
% acqmap = reshape(o{run}.sail{c}.acqMap(end).genes,625,10);
% t_points = out_of_sample_est(acqmap, optima, o{run}.reducedParSpace);
% scatter(t_points(:,1),t_points(:,2),'filled');
% legend('500PE Optima', 'cluster', '1000PE AcqMap');
% title('Acquisition map, linearly approx. reduced locations');
% 
% fig(6) = figure(6);hold off;
% scatter(o{run}.reducedParSpace(:,1),o{run}.reducedParSpace(:,2),'filled');
% hold on;
% scatter(memberLocations(:,1),memberLocations(:,2),'filled');
% predmap = reshape(o{run}.sail{c}.predMap(end).genes,625,10);
% t_points = out_of_sample_est(predmap, optima, o{run}.reducedParSpace);
% scatter(t_points(:,1),t_points(:,2),'filled');
% legend('500PE Optima', 'cluster', '1000PE PredMap');
% title('Prediction map, linearly approx. reduced locations');

save_figures(fig, './', ['sail_on_cluster_run' int2str(run) '_cluster' int2str(c) '_info_'], 12, [8 8]);


