function config = poemParamSet(mapDefaults,AEDefaults)
%POEMPARAMSET 

config.map                           = mapDefaults;
config.model                         = AEDefaults;
config.map.numInitSamples            = 256;

config.mutSelection                  = 0.1;
config.retryInvalid                  = true;

% Visualization and data management
config.display.illu                 = false;
config.display.illuMod              = 1;
end

