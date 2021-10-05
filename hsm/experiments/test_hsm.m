%% Train MLP
clear;clc;set(0,'DefaultFigureWindowStyle','docked');
tic;
%% Load testing dataset
%[x,t] = cho_dataset;t = t(1,:);
%load('notallvelos');x=params';t=cD';
%[x,t] = simplefit_dataset;
[x,t] = engine_dataset; t = t(1,:);

x = mapminmax(x); t = mapminmax(t);

%% Data configuration and normalization
data_cfg.num_inputs = size(x,1); data_cfg.num_outputs = size(t,1);
train_size = 0.3;
t = t'; x  = x'; cfg.num_outputs = size(t,2); [x,data_cfg.psx] = mapminmax(x'); [t,data_cfg.pst] = mapminmax(t');
numExp = 1;
    
%% Train HSM
train_errors = []; validation_errors = []; test_errors = [];
for i = 1:numExp
    dicethrow = rand(size(x,2),1);
    train_set = dicethrow < train_size; test_set = dicethrow >= train_size;
    nw{i}.cfg = params_hsm(x(:,train_set));
    nw{i} = train_hsm(x(:,train_set), t(:,train_set), nw{i});
end

%% Visualization
i = 1;
for j=1:50
    [prediction,layer_prediction,prediction_ids] = predict_hsm(x(:,test_set), nw{i});
end
toc

vis_sample_layers(nw{i},layer_prediction,prediction_ids,x(:,test_set),t(:,test_set));


