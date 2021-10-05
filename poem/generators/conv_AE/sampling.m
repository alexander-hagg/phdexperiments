function zSampled = sampling(encoderNet, x)
compressed = forward(encoderNet, x);
zMean = compressed;

sz = size(zMean);
zMean = reshape(zMean, [1,1,sz]);
zSampled = dlarray(zMean, 'SSCB');
end