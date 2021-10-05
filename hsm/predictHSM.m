%% predict_hsm - Get hierarchical surrogate model's activation.

function [prediction, layer_prediction, prediction_ids] = predictHSM( inputs, coordinates, model )
layers = get_layers(model);
layer_prediction = cell(model.cfg.hierarchy.num_layers,1);
layer_prediction{1} = predict_expert(inputs, model.hierarchy.Node{1}{1}); %% Get prediction of root model
prediction_ids{1} = ones(1,length(layer_prediction{1}));

% For every layer
for layer=2:model.cfg.hierarchy.num_layers
    layer_prediction{layer} = layer_prediction{layer-1};
    
    % Get all centers
    children = [];
    for parent=1:length(layers{layer-1})
        children = [children model.hierarchy.getchildren(layers{layer-1}(parent))];
    end
    if ~isempty(children)
        childModels = [model.hierarchy.Node{children}];
        childModels = reshape(childModels,3,length(childModels)/3);
        centers = [childModels{3,:}]';
        
        % Determine segment memberships
        if length(centers)==1
            idx = ones(1,size(inputs,2));
        elseif ~isempty(centers)
            distances = pdist2(coordinates,centers);
            if length(children) == 1
                idx = ones(1,size(inputs,2));
            else
                [~,idx] = sort(distances'); idx = idx(1,:);
            end
        end
        
        % Get residual activations from child models
        if ~isempty(centers)
            childIDs = unique(idx);
            for c = childIDs
                layer_prediction{layer}(:,idx==c) = layer_prediction{layer-1}(:,idx==c) - predict_expert(inputs(:,idx==c), childModels{1,children(c)-1});
            end
            prediction_ids{layer} = children(idx);
        end
    end
end

prediction = layer_prediction{end};