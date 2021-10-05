function output = getGroundTruth(fname, dname, varargin)
%% Get true results of each prediction map
%
% Syntax:  getGroundTruth(fname, dname, dummy)
%
% This function takes a domain folder and mat file name which contains:
%
%   output[1 X nRuns]
%       .predMap [1 X nSamplingPts]
%           .fitness
%           .cD
%           .cL
%       .d
%           .preciseEvaluate
%
% It then saves a file to disk with the same name with the suffix '.true'
% and the output struct to include
%
%   output[1 X nRuns]
%       .predMap [1 X nSamples]
%           .fitness_true
%           .cD_true
%           .cL_true
%
% Example:
%   fname = '~/Code/data/ffdSail/sailFFD';
%   dname = '~/Code/ffdSail/domains/foilFFD/';
%   getGroundTruth(fname,dname,true); % Test output
%
% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (BRSU)
% email: adam.gaier@h-brs.de
% Aug 2017; Last revision: 24-Aug-2017

%------------- BEGIN CODE --------------
% Read in input
parse = inputParser;
parse.addRequired('fname');
parse.addRequired('dname');
parse.addOptional('numSamples',[]);
parse.addOptional('dummy',false);
parse.parse(fname, dname, varargin{:});
iItr = parse.Results.numSamples; 
dummy = parse.Results.dummy;


addpath(genpath(dname)); load(fname);

% Determine number of maps to evaluate
if isempty(iItr); iItr = p.data.mapEvalSteps;end
nItr= size(iItr,2);

predMap = output.predMap;
tic;
for i = 1:nItr
    disp(['Iteration ' int2str(i) '/' int2str(nItr)]);
    genes = reshape(predMap(iItr(i)).genes,...
        [size(predMap(iItr(i)).genes,1) * ...
        size(predMap(iItr(i)).genes,2), d.dof]) ;
    
    % Just test if it's working?
    if dummy
        trueFit = predMap(iItr(i)).fitness(:)./10;
        trueVal = [predMap(iItr(i)).cD(:) predMap(iItr(i)).cL(:)]./10;
    else
        [trueVal, trueFit] = feval(d.preciseEvaluate,genes,d);
    end
    
    % Initialize the maps in the correct shape
    predMap(iItr(i)).fitness_true    = nan(size(predMap(iItr(i)).fitness));
    predMap(iItr(i)).cD_true         = nan(size(predMap(iItr(i)).cD));
    predMap(iItr(i)).cL_true         = nan(size(predMap(iItr(i)).cL));
    
    % Insert Values
    predMap(iItr(i)).fitness_true(:) = trueFit;
    predMap(iItr(i)).cD_true(:)      = trueVal(:,1);
    predMap(iItr(i)).cL_true(:)      = trueVal(:,2);
    
    disp(num2str(toc));
end
output.predMap = predMap;


output.samples.genes    = output.modelPred{1}.trainInput;
acqFunction = feval(d.createAcqFunction, output.modelPred, d);
[output.samples.fitness, predVals] = acqFunction(output.samples.genes);
output.samples.cD       = predVals{1}(:,1);
output.samples.cL       = predVals{2}(:,1);

shape = d.express(output.samples.genes);
area = squeeze(polyarea(shape(1,:,:), shape(2,:,:)));
areaPenalty = (1-(abs(area-d.base.area)./d.base.area)).^7;    
liftPenalty = ones(size(output.samples.genes,1),1);
liftPenalty(output.samples.cL < d.base.lift) = 1 / ((1-(output.samples.cL(output.samples.cL < d.base.lift)-d.base.lift)./d.base.lift).^2);

if dummy
    output.samples.cD_true      = output.samples.cD/10;
    output.samples.cL_true      = output.samples.cL/10;
    output.samples.fitness_true = output.samples.fitness/10;
else
    [trueVal, trueFit] = feval(d.preciseEvaluate,output.samples.genes,d);
    output.samples.cD_true      = trueVal(:,1);
    output.samples.cL_true      = trueVal(:,2);
    output.samples.fitness_true = trueFit;
end

%save([fname '_true'],'output','p','d')


disp(['Run done.'])

%------------- END OF CODE --------------






