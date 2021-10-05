%% get_segments - Get segments for a sample set

function [segment_ids, centers] = getSegments(cfg, coordinates)

[segment_ids,centers] = kmedoids(coordinates',cfg.kcenters);
segment_ids = segment_ids';
centers = centers';

end