%% train_hsm - Train a residual-hierarchical surrogate model
% Needs a model, preinitialized with @params_hsm, as well as samples
%
% A tree node in the hierarchy consists of a model, a logical array
% indicating the sample subset and the segment center

function [model] = trainHSM(inputs, targets, coordinates, model)
% Start hierarchy with a global model
selected = logical(ones(1,size(inputs,2)));
ROOT.cfg = model.cfg; ROOT.cfg.underfit = true;
hierarchy = tree({  train_expert(inputs(:,selected), targets(selected), ROOT.cfg),  ...
    selected, NaN                 });
layer_prediction{1} = hierarchy.Node{1}{1}.stats.train_prediction;

parents = [1]; % Set root node as parent

for layer=2:model.cfg.hierarchy.num_layers
    layer_prediction{layer} = layer_prediction{layer-1};
    newparents = [];
    for p=1:length(parents)
        samples = hierarchy.Node{parents(p)}{2};
        if sum(samples) > model.cfg.segmentation.min_samples
            [ids,centers] = getSegments(model.cfg.segmentation, coordinates(:,samples));
            for m=unique(ids)
                % Create logical array with segment's samples, but with original length of samples array
                selected = samples; selected(samples) = (ids==m);
                if sum(selected) > model.cfg.segmentation.min_samples;
                    parent = hierarchy.get(parents(p));
                    target = layer_prediction{layer-1}(selected)  - targets(selected); % Get residual from parent
                    [hierarchy node] = hierarchy.addnode(parents(p), ...
                        {train_expert(inputs(:,selected), target, model.cfg), ...
                        selected, centers(:,m)});
                    layer_prediction{layer}(selected) = layer_prediction{layer-1}(selected) - hierarchy.Node{node}{1}.stats.train_prediction;
                    newparents = [newparents node];
                end
            end
        end
    end
    parents = newparents;
end

model.hierarchy = hierarchy;
model.layer_prediction = layer_prediction;

end