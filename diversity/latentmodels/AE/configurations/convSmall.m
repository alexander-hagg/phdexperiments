function ae = convSmall(latentDim,resolution,numFilters)
%CONVSMALL Convolutional AE, small structure
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

