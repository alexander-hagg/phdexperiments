function reconstructionError = getReconstructionError(genomes,latent,manifold,getPhenotype)
%GETRECONSTRUCTIONERROR Summary of this function goes here
%   Detailed explanation goes here
latX(1,1,:,:) = cat(4,latent');
latX = dlarray(single(latX), 'SSCB');

xPred = [];
if size(latX,4) > 256
    subd = [1:256:size(latX,4) size(latX,4)];
    for i=1:length(subd)-1
        xPred = cat(4,xPred,sigmoid(forward(manifold.decoderNet, latX(:,:,:,subd(i):subd(i+1)))));
    end
else
    xPred = sigmoid(forward(manifold.decoderNet, latX));
end

[~,newPhenotypes] = getPhenotype(genomes);
xTrue = cat(4,newPhenotypes{:});
reconstructionError = zeros(size(genomes,1),1);
for i=1:size(xTrue,4)
    reconstructionError(i) = gather(extractdata(loss(xTrue(:,:,:,i), xPred(:,:,:,i), [], [])));
end
end

