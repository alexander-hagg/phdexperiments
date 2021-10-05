function [model,resultTupels] = trainGP(input,output,model,varargin)
%trainGP - Trains Gaussian Process model
% Given training input and output, optimizes given hyperParameters
%
% Syntax:  [d] = trainGP(input,output,d)
%
% Inputs:
%    input  - [samples X input dims]
%    output - [samples X 1]
%    d      - GP parameter struct
%
% Outputs:
%    GP_model - trained GP model
%
% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% May 2016; Last revision: 02-Aug-2016

%------------- BEGIN CODE --------------
gridsearch = true;
if nargin > 3
    gridsearch = varargin{1};
end

% Remove non-unique samples
[input,uniqueIDs] = unique(input,'rows','stable');
output = output(uniqueIDs);

model.hyp.mean = mean(output);

if gridsearch
    disp('Running grid search on length scale and signal variance');
    lengthscales = [-1 0 1];
    noiseVar = [-1 0 1];
    l = 1;
    for i=1:length(lengthscales)
        for j=1:length(noiseVar)
            hyp{l} = model.hyp;
            hyp{l}.cov(1) = lengthscales(i);
            hyp{l}.cov(2) = noiseVar(j);
            [hyp{l}, loglikelihoodEvolution] = minimize(hyp{l}, @gp, -model.functionEvals, @infExact, model.meanfunc, model.covfunc, model.likfunc, input, output);
            loglikelihood(l) = loglikelihoodEvolution(end);
            resultTupels(l,:) = [lengthscales(i), noiseVar(j), loglikelihoodEvolution(end)];
            l = l + 1;
        end
    end
    
    disp('Picking smallest loglikelihood');
    [minloglikelihood, IDloglikelihood] = min(loglikelihood);
    
    model.hyp = hyp{IDloglikelihood};
    model.loglikelihoodEvolution = minloglikelihood;
else
    disp('NOT running grid search on length scale and noise variance');
    [model.hyp, loglikelihoodEvolution] = minimize(model.hyp, @gp, -model.functionEvals, @infExact, model.meanfunc, model.covfunc, model.likfunc, input, output);
    model.loglikelihoodEvolution = loglikelihoodEvolution;
end


model.trainInput = input; model.trainOutput = output;

disp(['GP Training - Length scale found: ' num2str(model.hyp.cov(1))]);
disp(['GP Training - Noise var found: ' num2str(model.hyp.cov(2))]);

%------------- END OF CODE --------------
