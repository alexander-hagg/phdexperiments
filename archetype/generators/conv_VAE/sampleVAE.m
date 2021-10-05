function returnImage = sampleVAE(latentVector,decoderNet)
input = reshape(latentVector,1,1,size(latentVector,1),[]);
input = dlarray(input,'SSCB');
generatedImage = sigmoid(predict(decoderNet, input));
returnImage = gather(extractdata(generatedImage));
end
