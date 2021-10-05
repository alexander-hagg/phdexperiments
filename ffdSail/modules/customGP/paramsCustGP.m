function p = paramsCustGP(nInputs)
%paramsGP - Creates Gaussian Process parameter struct

%------------- BEGIN CODE --------------
p.name          = 'CustGP';

p.covfunc       = @covSEard;             
p.hyp.cov       = [zeros(nInputs,1);0]; % (unit vector in log space)

p.hyp.meandeg   = 2;
p.meanfunc      = {@meanPoly,p.hyp.meandeg};
p.hyp.mean      = [ones(p.hyp.meandeg*nInputs,1)];

p.likfunc       = @likGauss;     
p.hyp.lik       = log(0.1);

p.functionEvals = 100;      % function evals to optimize hyperparams


%------------- END OF CODE --------------