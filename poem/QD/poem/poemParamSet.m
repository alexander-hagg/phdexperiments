function p = poemParamSet(workDir)
%SAILPARAMSET infill configuration for surrogate-assistance

p.map                           = defaultParamSet;

p.retryInvalid                  = true;
p.numInitSamples                = 200;
p.numTotalSamples               = 2^9;
p.numIterations                 = 3;

%% Selection methods
% 'all'         - replace all manifold training samples with new map
% 'maxerror'       - add p.selectPerc % of most badly reconstructed solutions
% 'minerror'       - add p.selectPerc % of most badly reconstructed solutions
% 'random'      - add p.selectPerc % randomly from map
% 'combined'    - add p.selectPerc of most novel AND p.selectPerc of most
%                 badly reconstructed solutions

p.selectionMethod               = 'random'; 
p.selectPerc                    = 10;
%p.replacePerc                   = 10;

% Visualization and data management
p.display.illu              = false;
p.display.illuMod           = 25;
end

