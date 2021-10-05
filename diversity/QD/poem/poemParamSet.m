function p = poemParamSet(mapDefaults,AEDefaults)
%POEMPARAMSET 

p.map                           = mapDefaults;
p.model                         = AEDefaults;

p.retryInvalid                  = true;
p.numInitSamples                = 200;
p.numIterations                 = 8;

%% Selection methods
% 'all'         - replace all manifold training samples with new map
% 'maxerror'       - add p.selectPerc % of most badly reconstructed solutions
% 'minerror'       - add p.selectPerc % of most badly reconstructed solutions
% 'random'      - add p.selectPerc % randomly from map
% 'combined'    - add p.selectPerc of most novel AND p.selectPerc of most
%                 badly reconstructed solutions

p.selectionMethod               = 'all'; 
p.selectPerc                    = 100;
%p.replacePerc                   = 10;

% Visualization and data management
p.display.illu              = false;
p.display.illuMod           = 5;
end

