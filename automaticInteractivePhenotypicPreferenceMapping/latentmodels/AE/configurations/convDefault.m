function ae = convDefault(latentDim,resolution,numFilters)
%CONVDEFAULT Convolutional AE, default structure
%
imageSize = [resolution resolution 1];
weightsInitializer = 'glorot'; %he narrow-normal glorot
filterSize = 3;
stride = 2;

ae.encoderLG = layerGraph([
    imageInputLayer(imageSize,'Name','input_encoder','Normalization','none')
    convolution2dLayer(filterSize, numFilters, 'Padding','same', 'Stride', stride, 'Name', 'conv1', 'WeightsInitializer', weightsInitializer)
    %reluLayer('Name','relu1')
    convolution2dLayer(filterSize, numFilters, 'Padding','same', 'Stride', stride, 'Name', 'conv2', 'WeightsInitializer', weightsInitializer)
    %reluLayer('Name','relu2')
    fullyConnectedLayer(latentDim, 'Name', 'fc_encoder', 'WeightsInitializer', weightsInitializer)
    ]);

ae.decoderLG = layerGraph([
    imageInputLayer([1 1 latentDim],'Name','i','Normalization','none')
    transposedConv2dLayer(resolution/(stride*stride), numFilters, 'Cropping', 'same', 'Stride', resolution/(stride*stride), 'Name', 'transpose1', 'WeightsInitializer', weightsInitializer)
    %reluLayer('Name','relu1')
    transposedConv2dLayer(filterSize, numFilters, 'Cropping', 'same', 'Stride', stride, 'Name', 'transpose2', 'WeightsInitializer', weightsInitializer)
    %reluLayer('Name','relu2')
    transposedConv2dLayer(filterSize, numFilters, 'Cropping', 'same', 'Stride', stride, 'Name', 'transpose3', 'WeightsInitializer', weightsInitializer)
    %reluLayer('Name','relu3')
    transposedConv2dLayer(filterSize, 1, 'Cropping', 'same', 'Name', 'transpose4', 'WeightsInitializer', weightsInitializer)
    ]);

end

