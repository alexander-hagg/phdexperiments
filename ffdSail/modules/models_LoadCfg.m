function [modelParsAcq, modelParsPre, segmentMap, modelCfgs] = models_LoadCfg( d )
%MODELS_LOAD Summary of this function goes here
%   Detailed explanation goes here
d.alignedMap = true;
segmentMap = {d.express, d.categorize, d.featureMin, d.featureMax, d.alignedMap};

modelCfgs{1} = {'GP', 'GP', d.dof};
modelCfgs{2} = {'BN', 'Bootstrap', 'basemodel', 'MLP', 'useHierarchy', false, 'numMembers', 30};
modelCfgs{3} = {'BHMap', 'Bootstrap', 'basemodel', 'MLP', 'useHierarchy', true, 'numSegments', 4, 'numMembers', 30, 'mapPars', segmentMap};
modelCfgs{4} = {'BHPar', 'Bootstrap', 'basemodel', 'MLP', 'useHierarchy', true, 'numSegments', 4, 'numMembers', 30};

modelParsAcq = modelCfgs{1};
modelParsPre = modelCfgs{1};

end

