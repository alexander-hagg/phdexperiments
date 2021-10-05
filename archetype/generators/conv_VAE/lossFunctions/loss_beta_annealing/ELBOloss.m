function [elbo,reconstructionLoss,KL] = ELBOloss(x, xPred, zMean, zLogvar, epoch, numEpochs)
beta = 5;
squares = 0.5*(xPred-x).^2;
reconstructionLoss  = sum(squares, [1,2,3]);
KL = -.5 * sum(1 + zLogvar - zMean.^2 - exp(zLogvar), 1);

annealingScalar = linear_annealing(0, 1, epoch, numEpochs);
elbo = mean(reconstructionLoss + (annealingScalar+1) * beta * KL);
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

