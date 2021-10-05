function reconstructionLoss = loss(x, xPred, epoch, numEpochs)
squares = 0.5*(xPred-x).^2;
reconstructionLoss  = sum(squares, [1,2,3]);
reconstructionLoss = mean(reconstructionLoss);
