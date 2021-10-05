function p = poemParamSet(mapDefaults,AEDefaults)
%POEMPARAMSET 

p.map                           = mapDefaults;
p.model                         = AEDefaults;

p.retryInvalid                  = true;
p.numInitSamples                = 32;
p.numIterations                 = 2;

%% Selection methods
% 'all'         - replace all manifold training samples with new map
% 'maxerror'       - add p.selectPerc % of most badly reconstructed solutions
% 'minerror'       - add p.selectPerc % of most badly reconstructed solutions
% 'random'      - add p.selectPerc % randomly from map
% 'combined'    - add p.selectPerc of most novel AND p.selectPerc of most
%                 badly reconstructed solutions

p.selectionMethod               = 'all'; % 'all' 'ascend' 'descend'
p.selectionPerc                 = 0.75;

% Visualization and data management
p.display.illu              = true;
p.display.illuMod           = 1;
end

