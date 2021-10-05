function expert = train_expert(inputs, targets, cfg)

data.num_inputs = size(inputs,1); data.num_outputs = size(targets,1);

if ~isnan(cfg.pca_min_variance)
    [inputs,params_pca] = processpca(inputs,cfg.pca_min_variance);
end

if strcmp('MLP',cfg.experts.modelName)
    expert.cfg = feval(['params' cfg.experts.modelName],inputs, data);
    if isfield(cfg, 'underfit')
        expert.cfg.samples_per_dof = expert.cfg.samples_per_dof_forced_underfit; %UNDERFITTING
    end
elseif strcmp('GP',cfg.experts.modelName)
    expert.cfg = feval(['params' cfg.experts.modelName],size(inputs,2));
else
    disp(['ERROR: base model ' cfg.experts.modelName ' not implemented']);
end

if ~isnan(cfg.pca_min_variance)
    expert.cfg.pca = params_pca;
else
    expert.cfg.pca = [];    
end

if strcmp('MLP',cfg.experts.modelName)
    expert = trainMLP(inputs, targets, expert);
elseif strcmp('GP',cfg.experts.modelName)
    expert = trainGP(inputs, targets, expert);
end

expert.cfg.modelName = cfg.experts.modelName;
expert.cfg.data = data;

end