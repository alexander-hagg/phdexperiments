%runSail - Example usage script of sail function
% Running sail without arguments will return a hyperparameter struct of
% default values. These defaults can be changed in
% /sail/defaultParamSet.m
% 
% Running sail with a parameter struct as input will run the algorithm
%
%
% Other m-files required: defaultParamSet, sail, mapElites
% Other submodules required: gpml-wrapper
% 
%
% See also: mapElites, sail

% Author: Adam Gaier
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: adam.gaier@h-brs.de
% Nov 2016; Last revision: 01-Aag-2017

%------------- BEGIN CODE --------------
clear; 
nRuns = 1;
for iRep = {'parsec','ffd','cppn'}
    addpath(genpath('domains'));
    d = feval([iRep{:} '_Domain']);
    for iRun = 1:nRuns
        p = sail; % Get Default parameters
        
        % Adjust hyperparameters
        p.nTotalSamples     = 60;
        
        p.nGens = 250;
        
        p.data.mapEval      = true;
        p.data.mapEvalMod   = 50;
        % % % % % % % % % % % % % % %
        
        tstart = tic;        
        output{iRun} = sail(p,d);     
        runTime = toc(tstart);
        disp(['Run ' int2str(iRun) ' of ' int2str(nRuns)  ' completed in ' seconds2human(runTime)])
        disp(['Time remaining: ' seconds2human(runTime*(nRuns-iRun))])
        save(['sail_' iRep{:} '.mat']);
    end    
end









