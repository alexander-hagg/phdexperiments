function model = trainLatentModel(model,data,numEpochs,maxBatchSize,learnRate)
%TRAINVAE Summary of this function goes here
%   Detailed explanation goes here

encoderNet = dlnetwork(model.encoderLG);
decoderNet = dlnetwork(model.decoderLG);

executionEnvironment = "auto"; %auto cpu gpu

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
            adamupdate(decoderNet.Learnables, genGrad, avgGradientsDecoder, avgGradientsSquaredDecoder, iteration, learnRate);
        [encoderNet.Learnables, avgGradientsEncoder, avgGradientsSquaredEncoder] = ...
            adamupdate(encoderNet.Learnables, infGrad, avgGradientsEncoder, avgGradientsSquaredEncoder, iteration, learnRate);
    end
    elapsedTime = toc;
    
    % If we are testing at all
    if ~strcmp(data.testStore,'emptyStore')
        testDataTable = data.testStore.readall;
        XTest = cat(4,testDataTable.input{:,1});
        XTest = dlarray(single(XTest), 'SSCB');
        if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
            XTest = gpuArray(XTest);
        end
        
        [z, zMean, zLogvar] = sampling(encoderNet, XTest);
        xPred = sigmoid(forward(decoderNet, z));
        [loss,reconstructionLoss,regTerm] = ELBOloss(XTest, xPred, zMean, zLogvar, z, epoch, numEpochs);
        if epoch==1 || mod(epoch,50)==0
            disp("Epoch : " + epoch + ...
                " Test loss = "+gather(extractdata(loss))+ ...
                " Test reconstruction error = "+gather(extractdata(reconstructionLoss))+ ...
                " Test regularization term = "+gather(extractdata(regTerm))+ ...
                ". Time taken for epoch = "+ elapsedTime + "s")
        end
        model.losses(epoch) = gather(extractdata(loss));
    else
        if epoch==1 || mod(epoch,50)==0
            disp("Epoch : " + epoch + ...
                ". Time taken for epoch = "+ elapsedTime + "s");
            it=1;
            data.trainStore.reset;
            while hasdata(data.trainStore)
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
            
            
                [z, zMean, zLogvar] = sampling(encoderNet, XBatch);
                xPred = sigmoid(forward(decoderNet, z));
                [loss(it),reconstructionLoss(it),regTerm(it)] = ELBOloss(XBatch, xPred, zMean, zLogvar, z, epoch, numEpochs);
                it = it + 1;
            end
            disp("Epoch : " + epoch + ...
                " Mean Training loss = "+mean(gather(extractdata(loss)))+ ...
                " Mean Training reconstruction loss = "+mean(gather(extractdata(reconstructionLoss)))+ ...
                " Mean Training regularization term = "+mean(gather(extractdata(regTerm)))+ ...
                ". Time taken for epoch = "+ elapsedTime + "s")
        end
    end
end

model.encoderNet = encoderNet;
model.decoderNet = decoderNet;

% Normalize latent space
if ~strcmp(data.testStore,'emptyStore')
    allData = [readall(data.trainStore); readall(data.testStore)];
else
    allData = readall(data.trainStore);
end

model.latent = getLatent(allData{:,1},model);

if size(model.latent,1) > 2
    [coeff, score, latent, tsquared, explained] = pca(model.latent');
    model.latent = score(:,1:2);
    model.pcaCoeff = coeff(:,1:2);
else
    model.latent = model.latent';
end

[model.latent,model.normalization] = mapminmax(model.latent',0,1);
model.latent = model.latent';



end

