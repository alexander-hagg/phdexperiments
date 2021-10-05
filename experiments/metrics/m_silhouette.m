function val = m_silhouette(locations,clusterLabels)
%SILHOUETTE Summary of this function goes here
%   Detailed explanation goes here

eval = evalclusters(locations,clusterLabels+1,'silhouette');
val = eval.CriterionValues;
end

