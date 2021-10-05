function fig = viewPrototypes(optima, clusterIDs, estimatedLabels, d)
%VIEWPROTOTYPES Summary of this function goes here
%   Detailed explanation goes here

%% Overlapping shapes of every cluster
numPCExplained = 5; % Over how many principal components calculate "% of variance explained"
axisRanges = [0 1 -0.4 0.4];
for ii=1:length(clusterIDs)
    % Get all nD samples in this cluster
    samples = optima(estimatedLabels==clusterIDs(ii),:);
    for ss=1:size(samples,1)
        sampleExpressed{ii}(ss,:,:) = d.express(samples(ss,:));
    end
end

% Get medoids of clusters
for ii=1:length(clusterIDs)
    % Get all nD samples in this cluster
    samples = optima(estimatedLabels==clusterIDs(ii),:);
    [~,c] = kmedoids(samples,1);
    prototype{3}(ii,:,:) = d.express(c);
end

%
for ii=1:length(clusterIDs)
    % Get all nD samples in this cluster
    samples = optima(estimatedLabels==clusterIDs(ii),:);
    fig(ii) = figure; hold on;    
    for ss=1:size(samples,1)
        p1 = plot(squeeze(sampleExpressed{ii}(ss,1,:)),squeeze(sampleExpressed{ii}(ss,2,:)),'k','LineWidth',2);
        p1.Color(4) = min(max(5/size(samples,1),0),1);
    end
    p5 = plot(squeeze(prototype{3}(ii,1,:)),squeeze(prototype{3}(ii,2,:)),'g-','LineWidth',3);
    ax = gca;
    ax.XAxis.Visible =  'off';
    ax.YAxis.Visible =  'off';
    axis equal;axis(axisRanges);
    grid on;
end

end

