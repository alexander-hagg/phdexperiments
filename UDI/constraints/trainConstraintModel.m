function model = trainConstraintModel(coord,latent)
%TRAINCONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

% Configure GP model
%p.covfunc   = {@covMaternard, 3};
%p.hyp.cov   = [zeros(size(coord,2),1);0]; % (unit vector in log space)
%p.covfunc   = @covSEiso;
%p.hyp.cov   = [0;0]; % (unit vector in log space)
p.covfunc   = {@covMaterniso, 3};
p.hyp.cov   = [log(mean(mean(pdist2(coord,coord))));1]; % in log space: mean distance between samples

p.meanfunc  = {@meanConst};
p.likfunc   = @likGauss;
p.hyp.lik   = log(0.1);
p.functionEvals = 100;      % function evals to optimize hyperparams

latentMean = mean(latent,1);

% Train GP mapping
model.x = p;
model.x.hyp.mean = latentMean(1);
model.x.hyp = minimize_gpml(model.x.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, p.covfunc, p.likfunc, coord, latent(:,1));

model.y = p;
model.y.hyp.mean = latentMean(2);
model.y.hyp = minimize_gpml(model.y.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, p.covfunc, p.likfunc, coord, latent(:,2));

model.trainInput  = coord;
model.trainOutput = latent;


disp(['Length scales: ' num2str(model.x.hyp.cov(1)) ' - ' num2str(model.x.hyp.cov(2))]);
end

