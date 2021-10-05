function prediction = predictGP(gpModel, input)
%predictGP - Produces predictions of output of given inputs
% Given training input and output, parameters and hyperParameters, predicts
% output of new inputs, with mean and variance
%
% Syntax:  prediction = predictGP(gpModel, input)
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
trainInput  = gpModel.trainInput;
trainOutput = gpModel.trainOutput;

[m, s2] = gp(gpModel.hyp, @infExact, gpModel.meanfunc, gpModel.covfunc, gpModel.likfunc, ...
            trainInput, trainOutput, input);

% figure(2);hold off;
% scatter(1:length(m),m);
% hold on;
% scatter(1:length(s2),10*s2);
% scatter(1:length(s2),m+10*s2);
% 
% ax = gca;grid on;
% ax.YAxis.Limits = [0 1];
% legend('\mu', '\sigma', '\mu + \sigma');
% drawnow;

prediction(:,1) = m;
prediction(:,2) = s2;

%------------- END OF CODE --------------