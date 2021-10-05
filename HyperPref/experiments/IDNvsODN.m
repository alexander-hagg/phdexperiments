clear;
%% Configuration
addpath(genpath(pwd))                           % Set path to all modules
DOF = 16;DOMAIN = 'catmullRom';                 % Degrees of freedom, Catmull-Rom spline domain
d = domain(DOF);                                % Domain configuration
numLatentDims = 2;
m = cfgLatentModel('data/workdir',d.resolution, numLatentDims);% VAE configuration

% Create shapes variations
shapeParams = [1 0.3 1 0.3 1 0.3 1 0.3 zeros(1,8)];
numShapes = 20; scaling = [0.2 1.0]; rotation = [0 0.45*pi];
genomes = createShapeVariations(shapeParams,numShapes,scaling,rotation);

%% Three datasets: last two miss part of the training data
deSelect{1} = []; deSelect{2} = [8:13]; deSelect{3} = [14:20];
allGenomes{1} = reshape(genomes,[],d.dof);
allGenomes{2} = genomes;allGenomes{2}(deSelect{2},deSelect{2},:) = nan;allGenomes{2} = reshape(allGenomes{2},[],d.dof);allGenomes{2}(all(isnan(allGenomes{2})'),:) = [];
allGenomes{3} = genomes;allGenomes{3}(deSelect{3},deSelect{3},:) = nan;allGenomes{3} = reshape(allGenomes{3},[],d.dof);allGenomes{3}(all(isnan(allGenomes{3})'),:) = [];

% Produce phenotypes
x = 1:numShapes; y = x; [X,Y] = ndgrid(x,y);
for i=1:3
    % Get phenotypes
    [fitness{i},phenotypes{i}] = fitfun(allGenomes{i},d);
    
    % Adjust placement to missing shapes
    tX = X; tY = Y; tX(deSelect{i},deSelect{i}) = nan; tY(deSelect{i},deSelect{i}) = nan;
    placement = [tX(~isnan(tX)'),tY(~isnan(tY)')];    
    fig(i) = figure(i); hold off; ax = gca;
    showPhenotype(allGenomes{i},d,1.2,ax,placement); axis equal;
    axis([0 1.3*numShapes -1.3*numShapes 0]);
end

%% Train models
for i=1:3; allModels{i} = trainFeatures(phenotypes{i},m);end
save([DOMAIN '_bvae_1000_2.mat']);

%% Visualization
for i=1:3
    % Get predicted features for *all* shapes, incl. the ones missing from
    % training data for particular model.
    [~,phen] = getPhenotypeBoolean(phenotypes{i}, allModels{i}.encoderLG.Layers(1).InputSize(1));
    features{i} = getPrediction(phen,allModels{i});
    
    % Create latent samples
    minFeatures = 1.2*min(features{i}(:)); maxFeatures = 1.2*max(features{i}(:));
    nSamples = 10;
    x = minFeatures:(maxFeatures-minFeatures)/nSamples:maxFeatures; y = x; [X,Y] = ndgrid(x,y);
    varyCoords = [X(:),Y(:)]';
    input = []; input(1,1,:,:) = varyCoords; input = dlarray(input,'SSCB');
    genImgSample = sigmoid(predict(allModels{i}.decoderNet, input));
    genImgSample = gather(extractdata(genImgSample)); 
    
    % Get reproduced shapes from training data
    input = []; input(1,1,:,:) = features{i}'; input = dlarray(input,'SSCB');
    genImgTrain = sigmoid(predict(allModels{i}.decoderNet, input));
    genImgTrain = gather(extractdata(genImgTrain));    %reproduced
    % Place collected VAE outputs in latent space
    scale = d.resolution;
    [normVaryCoords,mapping] = mapminmax(varyCoords,-nSamples,nSamples);
    bitmapCoords = 1 + (ceil(scale*normVaryCoords)-min(ceil(scale*normVaryCoords(:))));
    imgSize = [0 -min(bitmapCoords(:)) + max(bitmapCoords(:))+scale];
    
    % Place missing training data, iff exists    
    if ~isempty(deSelect{i})
        missingGenomes = genomes(deSelect{i},deSelect{i},:); missingGenomes = reshape(missingGenomes,[],d.dof);
        [~,phen] = fitfun(missingGenomes,d);
        [~,phen] = getPhenotypeBoolean(phen, allModels{i}.encoderLG.Layers(1).InputSize(1));
        missingFeatures{i} = getPrediction(phen,allModels{i});
        input = []; input(1,1,:,:) = missingFeatures{i}'; input = dlarray(input,'SSCB');
        genImgMissing = sigmoid(predict(allModels{i}.decoderNet, input));
        genImgMissing = gather(extractdata(genImgMissing));        
    end
    
    % Place *original* shapes into latent space for visual confirmation
    [~,phen] = fitfun(allGenomes{1},d);
    [~,phen] = getPhenotypeBoolean(phen, allModels{i}.encoderLG.Layers(1).InputSize(1));
    allFeatures{i} = getPrediction(phen,allModels{i});
    
    
    
    %% Turn VAE outputs to viewable images
    % Create image with samples in latent coordinates
    clear img; img{1} = (zeros(range(imgSize),range(imgSize)));img{2} = img{1}; img{3} = img{1}; img{4} = img{1};
    for jj=1:size(genImgSample,4)
        coords = bitmapCoords(:,jj)';
        img{1}([coords(1):(coords(1)+scale-1)],[coords(2):(coords(2)+scale-1)]) = img{1}([coords(1):(coords(1)+scale-1)],[coords(2):(coords(2)+scale-1)]) + squeeze(genImgSample(:,:,:,jj));
    end
    img{1} = img{1} > 0.9;
    
    % Create image with training examples in latent coordinates
    normTrainCoords = mapminmax.apply(features{i}',mapping)';
    trainCoords = 1 + (ceil(scale*normTrainCoords)-min(ceil(scale*normVaryCoords(:))));
    for jj=1:size(genImgTrain,4)
        coords = trainCoords(jj,:);
        img{2}([coords(1):(coords(1)+scale-1)],[coords(2):(coords(2)+scale-1)]) = img{2}([coords(1):(coords(1)+scale-1)],[coords(2):(coords(2)+scale-1)]) + squeeze(genImgTrain(:,:,:,jj));
    end    
    img{2} = img{2} > 0.9;
    
    % Create image with missing training examples in latent coordinates
    if ~isempty(deSelect{i})
        normTrainCoords = mapminmax.apply(missingFeatures{i}',mapping)';
        trainCoords = 1 + (ceil(scale*normTrainCoords)-min(ceil(scale*normVaryCoords(:))));
        for jj=1:size(genImgMissing,4)
            coords = trainCoords(jj,:);
            img{3}([coords(1):(coords(1)+scale-1)],[coords(2):(coords(2)+scale-1)]) = img{3}([coords(1):(coords(1)+scale-1)],[coords(2):(coords(2)+scale-1)]) + squeeze(genImgMissing(:,:,:,jj));
        end
    end
    img{3} = img{3} > 0.9;
       
    % Create image with ground truth training examples in latent coordinates
    normTrainCoords = mapminmax.apply(allFeatures{i}',mapping)';
    trainCoords = 1 + (ceil(scale*normTrainCoords)-min(ceil(scale*normVaryCoords(:))));
    for jj=1:size(allFeatures{i},1)
        coords = trainCoords(jj,:);        
        img{4}([coords(1):(coords(1)+scale-1)],[coords(2):(coords(2)+scale-1)]) = img{4}([coords(1):(coords(1)+scale-1)],[coords(2):(coords(2)+scale-1)]) + phen{jj};
    end
    img{4} = img{4} > 0.9;
    
    % Add images together
    imgComplete = zeros(size(img{1}));
    imgComplete = 0.25*img{1};                              % Samples
    imgComplete(img{4}(:)>0) = 0.5*img{4}(img{4}(:)>0);         % Ground truth of all training examples
    imgComplete(img{2}(:)>0) = 0.75*img{2}(img{2}(:)>0);     % Selected training examples
    imgComplete(img{3}(:)>0) = img{3}(img{3}(:)>0);    % Missing training examples
    
    
    % Show image
    fig(99+i) = figure(99+i);
    hold off; ax=gca;
    imshow(imgComplete);
    hold on;
    colormap([1 1 1; 0 0 0; 0 0 1; 0 1 0; 1 0 0]);
    drawnow  
end



%%
save_figures(fig, '.', 'IDNODN', 12, [5 5])



