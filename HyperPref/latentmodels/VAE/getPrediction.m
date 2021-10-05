function [latent,xPred,xTrue] = getPrediction(phenotypes,model)
%GETPREDICTION Summary of this function goes here
%   Detailed explanation goes here
xTrue = cat(4,phenotypes);
latent = getLatent(phenotypes,model);
latX(1,1,:,:) = cat(4,latent);
latX = dlarray(single(latX), 'SSCB');
latent = latent';

if nargin == 1
    return;
end

result = [];
if size(latX,4) > 256
    subd = [0:256:size(latX,4) size(latX,4)];
    for i=1:length(subd)-1
        result = cat(4,result,sigmoid(forward(model.decoderNet, latX(:,:,:,subd(i)+1:subd(i+1)))));
    end
else
    result = sigmoid(forward(model.decoderNet, latX));
end
xPred = gather(extractdata(result));
end

