function reconstructionLoss = loss(x, xPred, epoch, numEpochs)
reconstructionLoss = crossentropy(xPred,x);
