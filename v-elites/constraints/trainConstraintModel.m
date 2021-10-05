function model = trainConstraintModel(samples,latent)
%trainConstraintModel - train constraint model
%
% Syntax:  model = trainConstraintModel(coord,latent)
%
% Inputs:
%    samples - coordinates in original space
%    latent - coordinates in latent space (t-SNE)
%
% Outputs:
%    model
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Aug 2019; Last revision: 15-Aug-2019
%
%------------- BEGIN CODE --------------

% Configure GP model
p.covfunc   = {@covMaterniso, 3};
p.hyp.cov   = [log(mean(mean(pdist2(samples,samples))));1]; % in log space: mean distance between samples

p.meanfunc  = {@meanConst};
p.likfunc   = @likGauss;
p.hyp.lik   = log(0.1);
p.functionEvals = 100;      % function evals to optimize hyperparams

latentMean = mean(latent,1);

% Train GP mapping
model.x = p;
model.x.hyp.mean = latentMean(1);
model.x.hyp = minimize_gpml(model.x.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, p.covfunc, p.likfunc, samples, latent(:,1));

model.y = p;
model.y.hyp.mean = latentMean(2);
model.y.hyp = minimize_gpml(model.y.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, p.covfunc, p.likfunc, samples, latent(:,2));

model.trainInput  = samples;
model.trainOutput = latent;

end

%------------- END CODE --------------
