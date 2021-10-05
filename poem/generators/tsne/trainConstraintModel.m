function m = trainConstraintModel(samples,m)
%trainConstraintm.model - train constraint m.model
%
% Syntax:  m.model = trainConstraintm.model(coord,latent)
%
% Inputs:
%    samples - coordinates in original space
%    latent - coordinates in latent space (t-SNE)
%
% Outputs:
%    m.model
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Aug 2019; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

% Run t-SNE
latent = fast_tsne(samples, m.numDims_DR, m.numDims_ORG, m.perplexity, m.speedQualitytradeOff);

% Configure GP m.model
p.covfunc   = {@covMaterniso, 3};
p.hyp.cov   = [log(mean(mean(pdist2(samples,samples)))); 1]; % in log space: mean distance between samples, noise

p.meanfunc  = {@meanConst};
p.likfunc   = @likGauss;
p.hyp.lik   = log(0.01);
p.functionEvals = 100;      % function evals to optimize hyperparams

latentMean = mean(latent,1);

% Train GP mapping
%%
m.model.x = p;
m.model.x.hyp.mean = latentMean(1);
[m.model.x.hyp, m.model.x.stats.fX,i] = minimize_gpml(m.model.x.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, p.covfunc, p.likfunc, samples, latent(:,1));

m.model.y = p;
m.model.y.hyp.mean = latentMean(2);
[m.model.y.hyp, m.model.y.stats.fX] = minimize_gpml(m.model.y.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, p.covfunc, p.likfunc, samples, latent(:,2));

m.model.trainInput  = samples;
m.model.trainOutput = latent;

end

%------------- END CODE --------------
