function model = trainConstraints(coord,latent)
%TRAINCONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

% Configure GP model
p.covfunc   = {@covMaternard, 3};
p.hyp.cov   = [zeros(size(coord,2),1);0]; % (unit vector in log space)
p.meanfunc  = {@meanConst};
p.hyp.mean  = 0;
p.likfunc   = @likGauss;
p.hyp.lik   = log(0.1);
p.functionEvals = 100;      % function evals to optimize hyperparams

% Train GP mapping
model.x.hyp = minimize_gpml(p.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, ...
    p.covfunc, p.likfunc, coord, latent(:,1));
model.y.hyp = minimize_gpml(p.hyp,@gp, -p.functionEvals, @infExact, p.meanfunc, ...
    p.covfunc, p.likfunc, coord, latent(:,2));

model.p = p;
model.trainInput  = coord;
model.trainOutput = latent;
end

