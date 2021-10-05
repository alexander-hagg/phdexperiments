function [elbo,reconstructionLoss,KL] = ELBOloss(x, xPred, zMean, zLogvar, z, epoch, numEpochs)
beta = 10;
maxNats = 1; % not realyl calculating nats
%squares = 0.5*(xPred-x).^2;
%reconstructionLoss  = sum(squares, [1,2,3]);
reconstructionLoss = crossentropy(xPred,x,'TargetCategories','independent');

KL = -.5 * sum(1 + zLogvar - zMean.^2 - exp(zLogvar), 1);

annealingScalar = linear_annealing(0, maxNats, epoch, numEpochs);
elbo = mean(reconstructionLoss + beta * (KL - annealingScalar));

% For output
reconstructionLoss = mean(reconstructionLoss);
KL = mean(KL);

function anneal = linear_annealing(init, fin, step, maxSteps)

if step == 0
    anneal = fin;
    return;
end

if fin <= init
    disp('Annealing: final value is not larger than initial value!');
    anneal = fin;
    return;
end

delta = fin - init;
anneal = init + step*(delta/maxSteps);

