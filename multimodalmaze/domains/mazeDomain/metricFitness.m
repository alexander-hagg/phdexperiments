function trajectoryDistances = metricFitness(trajectories)
%METRICFITNESS Summary of this function goes here
%   Detailed explanation goes here
    segments = abs(trajectories(:,2:end,1:2)-trajectories(:,1:end-1,1:2));
    trajectoryDistances = sum(sqrt(segments(:,:,1).^2 + segments(:,:,2).^2),2);
end

