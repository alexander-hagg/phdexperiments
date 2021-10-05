function [GP_model] = trainGP(input,output,d)
%trainGP - Trains Gaussian Process model
% Given training input and output, optimizes given hyperParameters
%
% Syntax:  [output1,output2] = function_name(input1,output,gpParams)
%
% Inputs:
%    input  - [samples X input dims]
%    output - [samples X 1]
%    d      - GP parameter struct
%

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% May 2016; Last revision: 02-Aug-2016

%------------- INPUT PARSING -----------
parse = inputParser;
parse.addRequired('input');
parse.addRequired('output');
parse.addRequired('d');
parse.parse(input,output,d);
%functionEvals = parse.Results.functionEvals;

% Remove non-unique samples
[input,uniqueIDs] = unique(input,'rows');
output = output(uniqueIDs);

%------------- BEGIN CODE --------------
GP_model.hyp = minimize_gpml(d.hyp,@gp, -d.functionEvals, @infExact, d.meanfunc, ...
               d.covfunc, d.likfunc, input, output);
GP_model.trainInput = input;
GP_model.trainOutput= output;
GP_model.params     = d;


%------------- END OF CODE --------------