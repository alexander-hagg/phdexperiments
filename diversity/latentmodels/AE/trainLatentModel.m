function model = trainLatentModel(model,data,numEpochs,maxBatchSize,learnRate)
%TRAINLATENTMODEL Summary of this function goes here
%   Detailed explanation goes here

encoderNet = dlnetwork(model.encoderLG);
decoderNet = dlnetwork(model.decoderLG);

executionEnvironment = "gpu"; %auto cpu gpu

numTrainImages = data.trainStore.NumObservations;

miniBatchSize = min(maxBatchSize,numTrainImages);
data.trainStore.MiniBatchSize = miniBatchSize;

iteration = 0;

avgGradientsEncoder = [];
avgGradientsSquaredEncoder = [];
avgGradientsDecoder = [];
avgGradientsSquaredDecoder = [];


for epoch = 1:numEpochs
    tic;
    % Reset and shuffle datastore.
    reset(data.trainStore);
    data.trainStore = shuffle(data.trainStore);
    while hasdata(data.trainStore)
        iteration = iteration + 1;
        
        % Read mini-batch of data.
        batchDataTable = read(data.trainStore);
        
        % Ignore last partial mini-batch of epoch.
        if size(batchDataTable,1) < miniBatchSize
            continue
        end
        
        XBatch = cat(4,batchDataTable.input{:,1});
        XBatch = dlarray(single(XBatch), 'SSCB');
        
        if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
            XBatch = gpuArray(XBatch);
        end
        
        [infGrad, genGrad] = dlfeval(...
            @modelGradients, encoderNet, decoderNet, XBatch, epoch, numEpochs);
        
        [decoderNet.Learnables, avgGradientsDecoder, avgGradientsSquaredDecoder] = ...
            adamupdate(decoderNet.Learnables, ...
            genGrad, avgGradientsDecoder, avgGradientsSquaredDecoder, iteration, learnRate);
        [encoderNet.Learnables, avgGradientsEncoder, avgGradientsSquaredEncoder] = ...
            adamupdate(encoderNet.Learnables, ...
            infGrad, avgGradientsEncoder, avgGradientsSquaredEncoder, iteration, learnRate);
    end
    elapsedTime = toc;
    testDataTable = data.testStore.readall;
    XTest = cat(4,testDataTable.input{:,1});
    XTest = dlarray(single(XTest), 'SSCB');
    if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
        XTest = gpuArray(XTest);
    end
    
    z = sampling(encoderNet, XTest);
    xPred = sigmoid(forward(decoderNet, z));
    reconstructionLoss = loss(XTest, xPred, epoch, numEpochs);
    
    disp("Epoch : " + epoch + ...
        " Test reconstruction loss = "+gather(extractdata(reconstructionLoss))+ ...
        ". Time taken for epoch = "+ elapsedTime + "s")
    model.losses(epoch) = gather(extractdata(reconstructionLoss));
end

model.encoderNet = encoderNet;
model.decoderNet = decoderNet;

% Normalize latent space
allData = [readall(data.trainStore); readall(data.testStore)];
model.latent = getLatent(allData{:,1},model);
[~,model.normalization] = mapminmax(model.latent,0,1);

end

