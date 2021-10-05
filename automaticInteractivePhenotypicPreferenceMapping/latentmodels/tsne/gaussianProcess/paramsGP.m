function p = paramsGP(nInputs)
%paramsGP - Creates Gaussian Process parameter struct
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    nInputs - DOF of representation, number of parameters
%
% Outputs:
%    p - struct - GP configuration
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Jun 2016; Last revision: 02-Jun-2016

%------------- BEGIN CODE --------------
p.name = 'GP';

p.covfunc   = @covSEiso;             
p.hyp.cov   = [1;0]; % (unit vector in log space)

p.meanfunc  = {@meanConst};  
p.hyp.mean  = 0;

p.likfunc   = @likGauss;     
p.hyp.lik   = log(0.1);

p.functionEvals = 100;      % function evals to optimize hyperparams


%------------- END OF CODE --------------