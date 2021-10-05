function [elbo,crossentropy,KL] = ELBOloss(x, xPred, zMean, zLogvar, epoch, numEpochs)
% Binary cross entropy for bernoulli distribution (treat image 0's and 1's as classification problem)

N = size(x,1)*size(x,2)*size(x,3);
%rX = reshape(x,[],size(x,4));
%rPX = reshape(xPred,[],size(xPred,4));
rX = x(:);
rPX = xPred(:);

innerTerm = rX.*log(rPX) + (1 - rX).*log(1 - rPX);
crossentropy = -(1/N) * sum(innerTerm);

KL = -.5 * sum(1 + zLogvar - zMean.^2 - exp(zLogvar), 1);
KL = sum(KL);

elbo = crossentropy + KL;

%beta = 5;
%annealingScalar = linear_annealing(0, 1, epoch, numEpochs);
%elbo = mean(crossentropy + annealingScalar * beta * KL);


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