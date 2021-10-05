function infill = infillParamSet(varargin)
%SAILPARAMSET infill configuration for surrogate-assistance
    
infill.nTotalSamples        = 100;
infill.nAdditionalSamples   = 10;
infill.retryInvalid         = false;
infill.intermediateMaps     = true;

infill.display.illu              = false;
infill.display.illuMod           = 256;

end

