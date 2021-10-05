modelNames = {'Bootstrap', 'Bootstrap', 'Bootstrap', 'GP'};
paramNames = {'BANN', 'BHSM_MLP_map', 'BHSM_MLP_par', 'GP'};
d = ffd_Domain;%d = velo_Domain; %d = af_Domain;
xp = 2;
d.surrogateName = modelNames{xp};
d.paramNames = paramNames;
d.modelNum = xp;

% Set parameters, including specific for HSM.
for i=1:2 % lift and drag
    d.params{i}= feval(['params' d.surrogateName], d.dof);
    d.params{i}.surrogateName = d.surrogateName;
    if strcmp(paramNames{xp},'BHSM_MLP_par')||strcmp(paramNames{xp},'BHSM_MLP_map')
        d.params{i}.memberstr = 'HSM'; d.params{i}.hsmmember = 'MLP';
    end
    d.params{i}.express = d.express;
    d.params{i}.categorize = d.categorize;
    d.params{i}.featureMin = d.featureMin;
    d.params{i}.featureMax = d.featureMax;
    if strcmp(paramNames{xp},'BHSM_MLP_map');d.params{i}.useFeatureMap = true;end
    if strcmp(paramNames{xp},'BANN')||strcmp(paramNames{xp},'BHSM_MLP_par')||strcmp(paramNames{xp},'BHSM_MLP_map')
        d.varCoef = 1;% Adjust variance weight to deal with order of magnitude variance discrepancy between GP and bootstrapped models
    end
end
%%
in = rand(1000,10);
intest = randn(100000,10);
out = sin(5*in(:,1)) + 0.01*randn(1000,1);
%%
tic
%for i=1:1
%model = trainBootstrap(in,out,d.params{1});
%end
pred = predictBootstrap(model,intest);
toc
