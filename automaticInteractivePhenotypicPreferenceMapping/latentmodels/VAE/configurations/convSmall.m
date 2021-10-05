function vae = convSmall(latentDim,resolution,numFilters)
%CONFIGUREVAE VAE is configured as in 
% Burgess, C. P., Higgins, I., Pal, A., Matthey, L., Watters, N., Desjardins, G., & Lerchner, A. (2017). 
% Understanding disentangling in ? -VAE
% 10 Apr 2018, (Nips).
%
imageSize = [resolution resolution 1];

vae.encoderLG = layerGraph([
    imageInputLayer(imageSize,'Name','input_encoder','Normalization','none')
    convolution2dLayer(3, numFilters, 'Padding','same', 'Stride', 2, 'Name', 'conv1')
    reluLayer('Name','relu1')    
    fullyConnectedLayer(2 * latentDim, 'Name', 'fc_encoder')
    ]);

vae.decoderLG = layerGraph([
    imageInputLayer([1 1 latentDim],'Name','i','Normalization','none')
    transposedConv2dLayer(resolution/4, numFilters, 'Cropping', 'same', 'Stride', resolution/4, 'Name', 'transpose1')
    reluLayer('Name','relu1')
    transposedConv2dLayer(3, numFilters, 'Cropping', 'same', 'Stride', 4, 'Name', 'transpose3')
    reluLayer('Name','relu3')
    transposedConv2dLayer(3, 1, 'Cropping', 'same', 'Name', 'transpose4')
    ]);

end

