%% predict_expert - Get prediction from expert model in the hierarchy

function prediction = predict_expert(inputs, model)

if isfield(model.cfg, 'pca') && ~isempty(model.cfg.pca)
    inputs = processpca('apply',inputs,model.cfg.pca);
end

if strcmp('MLP',model.cfg.modelName)
    prediction = predictMLP(inputs,model)';
elseif strcmp('GP',model.cfg.modelName)
    prediction = predictGP(inputs,model)';
else
    disp(['ERROR: base model ' model.cfg.modelName ' not implemented']);
end
end