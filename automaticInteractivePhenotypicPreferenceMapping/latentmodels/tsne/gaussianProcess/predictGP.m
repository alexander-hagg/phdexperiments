function prediction = predictGP(gpModel, input)
%predictGP - Produces predictions of output of given inputs
% Given training input and output, parameters and hyperParameters, predicts
% output of new inputs, with mean and variance
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    gpModel - trained GP model
%    input - locations for predictions
%
% Outputs:
%    prediction - prediction for locations 'input'
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
% May 2016; Last revision: 16-May-2016

%------------- BEGIN CODE --------------
hyp = gpModel.hyp;
trainInput  = gpModel.trainInput;
trainOutput = gpModel.trainOutput;
p = gpModel.params;

[m, s2] = gp(hyp, @infExact, p.meanfunc, p.covfunc, p.likfunc, ...
            trainInput, trainOutput, input);

prediction(:,1) = m;
prediction(:,2) = s2;

%------------- END OF CODE --------------