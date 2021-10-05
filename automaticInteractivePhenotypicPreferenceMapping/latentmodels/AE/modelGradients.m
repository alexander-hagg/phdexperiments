function [infGrad, genGrad] = modelGradients(encoderNet, decoderNet, x, epoch, numEpochs)
z = sampling(encoderNet, x);
xPred = sigmoid(forward(decoderNet, z));
losses = loss(x, xPred, epoch, numEpochs);
[genGrad, infGrad] = dlgradient(losses, decoderNet.Learnables, ...
    encoderNet.Learnables);
end