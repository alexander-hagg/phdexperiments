function zMean = getVAELatent(vae,XTest)
if iscell(XTest)
    XBatch = cat(4,XTest{:});
else
    XBatch = cat(4,XTest);
end
XBatch = dlarray(single(XBatch), 'SSCB');

%if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
%    XBatch = gpuArray(XBatch);
%end

[~, zMean] = sampling(vae.encoderNet, XBatch);
zMean = stripdims(zMean)';
zMean = gather(extractdata(zMean));
zMean = zMean';
end
