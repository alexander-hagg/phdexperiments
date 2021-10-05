function vae = configureVAE_default(latentDim,resolution,numFilters)
%CONFIGUREVAE VAE is configured as in 
% Burgess, C. P., Higgins, I., Pal, A., Matthey, L., Watters, N., Desjardins, G., & Lerchner, A. (2017). 
% Understanding disentangling in ? -VAE
% 10 Apr 2018, (Nips).
%
imageSize = [resolution resolution 1];

vae.encoderLG = layerGraph([
    imageInputLayer(imageSize,'Name','input_encoder','Normalization','none')
    convolution2dLayer(3, numFilters, 'Padding','same', 'Stride', 2, 'Name', 'conv1', 'WeightsInitializer', 'he')
    reluLayer('Name','relu1')
    convolution2dLayer(3, numFilters, 'Padding','same', 'Stride', 2, 'Name', 'conv2', 'WeightsInitializer', 'he')
    reluLayer('Name','relu2')
    fullyConnectedLayer(2 * latentDim, 'Name', 'fc_encoder', 'WeightsInitializer', 'he')
    ]);

vae.decoderLG = layerGraph([
    imageInputLayer([1 1 latentDim],'Name','i','Normalization','none')
    transposedConv2dLayer(resolution/4, numFilters, 'Cropping', 'same', 'Stride', resolution/4, 'Name', 'transpose1', 'WeightsInitializer', 'he')
    reluLayer('Name','relu1')
    transposedConv2dLayer(3, numFilters, 'Cropping', 'same', 'Stride', 2, 'Name', 'transpose2', 'WeightsInitializer', 'he')
    reluLayer('Name','relu2')
    transposedConv2dLayer(3, numFilters, 'Cropping', 'same', 'Stride', 2, 'Name', 'transpose3', 'WeightsInitializer', 'he')
    reluLayer('Name','relu3')
    transposedConv2dLayer(3, 1, 'Cropping', 'same', 'Name', 'transpose4', 'WeightsInitializer', 'he')
    ]);

end

