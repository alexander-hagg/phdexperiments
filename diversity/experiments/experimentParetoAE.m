% Ground Truth Pareto Set
d.resolution = 64;
nShapes = 20
PARETO_genomes = [];
axDev = d.ranges(1,1):(range(d.ranges(1,:))./(nShapes-1)):d.ranges(1,2);
radDev = d.ranges(9,1):(range(d.ranges(9,:))./(nShapes-1)):d.ranges(9,2);
for i=1:nShapes
    for j=1:nShapes
        PARETO_genomes(end+1,:) = [axDev(i).*ones(1,DOF/2), radDev(j).*ones(1,DOF/2)];
    end
end
[PARETO_fitness,PARETO_polygons,~] = d.fitfun(PARETO_genomes);
PARETO_features = categorize(PARETO_polygons, d)';
[~,PARETO_phenotypes] = getPhenotypeBoolean(PARETO_polygons,d.resolution);

%%
p = defaultParamSet;
poemCfg = poemParamSet(p,getAEConfig('data/workdir',d.resolution));
poemCfg.model = poemCfg.model.train(PARETO_phenotypes);

%%

[latent,xPred,xTrue] = getPrediction(PARETO_phenotypes,poemCfg.model);
figure(2);scatter(latent(:,1),latent(:,2));grid on;axis equal;

figure(1); hold off;
for ii=1:length(PARETO_phenotypes)
    subplot(1,2,1);
    imagesc(PARETO_phenotypes{ii});
    subplot(1,2,2);
    imagesc(squeeze(xPred(:,:,:,ii)));
    title(ii);
    drawnow;
    %pause(0.01);
end

%% Sample latent space
% Sample latent space equidistantly
minLat = min(latent);
maxLat = max(latent);

N = 20;
x = linspace(minLat(1), maxLat(1), N);
y = linspace(minLat(2), maxLat(2), N);

[X,Y] = ndgrid(x,y);

samplePts = [X(:),Y(:)];

returnImage = sampleAE(samplePts',poemCfg.model.decoderNet);
returnImage = squeeze(returnImage);
returnImage = permute(returnImage,[3,1,2]);

clear out;
for i=1:size(returnImage,1)
    out{i} = squeeze(returnImage(i,:,:));
end

figs(5) = figure(5);hold off;
out = imtile(out,'GridSize', [N N]);
imagesc(out);

%% Show shapes in learned latent space

latent = getPrediction(PARETO_phenotypes,poemCfg.model);
normLatent = mapminmax(latent',0,1); normLatent = 50*normLatent';
figs(3) = figure(3);hold off;
showPhenotype(PARETO_genomes,d,figs(3),normLatent);

%% Test feature code
poemCfg.model = trainAEfeatures(PARETO_phenotypes,poemCfg.model);
features = predictAEfeatures(PARETO_phenotypes,poemCfg.model)
features = 50*features;

scatter(features(:,1),features(:,2));




%%end code