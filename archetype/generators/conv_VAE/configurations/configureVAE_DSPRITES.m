function vae = configureVAE_DSPRITES(latentDim,resolution,numFilters)
%CONFIGUREVAE VAE is configured as in 
% Burgess, C. P., Higgins, I., Pal, A., Matthey, L., Watters, N., Desjardins, G., & Lerchner, A. (2017). 
% Understanding disentangling in ? -VAE
% 10 Apr 2018, (Nips).
%
imageSize = [resolution resolution 1];

vae.encoderLG = layerGraph([
    imageInputLayer(imageSize,'Name','input_encoder','Normalization','none')
    convolution2dLayer(4, numFilters, 'Padding','same', 'Stride', 2, 'Name', 'conv0')
    reluLayer('Name','relu0')
    convolution2dLayer(4, numFilters, 'Padding','same', 'Stride', 2, 'Name', 'conv1')
    reluLayer('Name','relu1')
    convolution2dLayer(4, numFilters, 'Padding','same', 'Stride', 2, 'Name', 'conv2')
    reluLayer('Name','relu2')
    convolution2dLayer(4, numFilters, 'Padding','same', 'Stride', 2, 'Name', 'conv3')
    reluLayer('Name','relu3')
    fullyConnectedLayer(256, 'Name', 'fc_encoder0')
    reluLayer('Name','relu4')
    fullyConnectedLayer(256, 'Name', 'fc_encoder1')
    reluLayer('Name','relu5')    
    fullyConnectedLayer(2 * latentDim, 'Name', 'fc_encoder')
    ]);

vae.decoderLG = layerGraph([
    imageInputLayer([1 1 latentDim],'Name','i','Normalization','none')  
    transposedConv2dLayer(resolution/8, numFilters, 'Cropping', 'same', 'Stride', resolution/8, 'Name', 'transpose0')
    reluLayer('Name','relu0')
    transposedConv2dLayer(4, numFilters, 'Cropping', 'same', 'Stride', 2, 'Name', 'transpose1')
    reluLayer('Name','relu1')
    transposedConv2dLayer(4, numFilters, 'Cropping', 'same', 'Stride', 2, 'Name', 'transpose2')
    reluLayer('Name','relu2')
    transposedConv2dLayer(4, numFilters, 'Cropping', 'same', 'Stride', 2, 'Name', 'transpose3')
    reluLayer('Name','relu3')
    transposedConv2dLayer(3, 1, 'Cropping', 'same', 'Name', 'transpose4')
    ]);

end

