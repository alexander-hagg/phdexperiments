function val = m_calinski(locations,clusterLabels)
%CALINSKI Summary of this function goes here
%   Detailed explanation goes here

eval = evalclusters(locations,clusterLabels+1,'silhouette');
val = eval.CriterionValues;
end

