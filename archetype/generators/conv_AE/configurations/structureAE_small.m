function ae = structureAE_small(latentDim,resolution,numFilters)
%CONFIGUREVAE AE is configured as in 
% Burgess, C. P., Higgins, I., Pal, A., Matthey, L., Watters, N., Desjardins, G., & Lerchner, A. (2017). 
% Understanding disentangling in ? -VAE
% 10 Apr 2018, (Nips).
%
imageSize = [resolution resolution 1];
filterSize = 3;
stride = 2;

ae.encoderLG = layerGraph([
    imageInputLayer(imageSize,'Name','input_encoder','Normalization','none')
    convolution2dLayer(filterSize, numFilters, 'Padding','same', 'Stride', stride, 'Name', 'conv1')
    %reluLayer('Name','relu1')    
    fullyConnectedLayer(latentDim, 'Name', 'fc_encoder')
    ]);

ae.decoderLG = layerGraph([
    imageInputLayer([1 1 latentDim],'Name','i','Normalization','none')
    transposedConv2dLayer(resolution/stride, numFilters, 'Cropping', 'same', 'Stride', resolution/stride, 'Name', 'transpose1')
    %reluLayer('Name','relu1')
    transposedConv2dLayer(filterSize, numFilters, 'Cropping', 'same', 'Stride', stride, 'Name', 'transpose3')
    %reluLayer('Name','relu3')
    transposedConv2dLayer(3, 1, 'Cropping', 'same', 'Name', 'transpose4')
    ]);

end

