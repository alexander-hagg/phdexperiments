function [ layers ] = get_layers( model )
%GET_LAYERS Summary of this function goes here
%   Detailed explanation goes here
layers{1} = [1];
for l=2:model.cfg.hierarchy.num_layers
    layers{l} = find(ismember(model.hierarchy.Parent,layers{l-1}));
end

end

