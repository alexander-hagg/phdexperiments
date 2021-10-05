function [output,reconstructionError] = getAELatent(XTest,model)
if iscell(XTest)
    XBatch = cat(4,XTest{:});
else
    XBatch = cat(4,XTest);
end
XBatch = dlarray(single(XBatch), 'SSCB');

%if (executionEnvironment == "auto" && canUseGPU) || executionEnvironment == "gpu"
%    XBatch = gpuArray(XBatch);
%end

% Divide up dataset to prevent out of memory on GPU
if length(XTest) > 1024
    subd = [1:1024:length(XTest) length(XTest)];
    for i=1:length(subd)-1
        zSampled = sampling(model.encoderNet, XBatch(:,:,:,subd(i):subd(i+1)));
        zSampled = squeeze(stripdims(zSampled))';
        zSampled = gather(extractdata(zSampled));
        output(subd(i):subd(i+1),:) = zSampled;
    end
else
    output = sampling(model.encoderNet, XBatch);
    output = squeeze(stripdims(output))';
    output = gather(extractdata(output));    
end
output = output';

%if isfield(model,'normalization')
%    output = mapminmax('apply',output,model.normalization);
%end

%output = double(output);
end
