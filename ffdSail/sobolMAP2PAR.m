systemInit;

d = ffd_Domain;

clrsmap = [1:5:25];
edges = {0:0.04:1 0:0.04:1};
sobSequence  = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');
observation = sobSequence(1:(1+1000)-1,:);

map = parula(35);   % a 35 element colormap
I = zeros(25);      % blank image
I(1) = 1;           % upper left corner is 1
I = bwdist(I);      % distance to upper left corner

clrsmap = repelem(clrsmap,1,5);
binColorV = repmat(clrsmap,25,1)/25;
binColorU = binColorV';
binColorW = ones(size(binColorV));
binColors = cat(3,binColorU,binColorV,binColorW);

fig = [];

fig(end+1) = figure(100);hold off;
im = imshow(binColors);
title('Bins in Feature Space');
%%

%dimReduxMethods = {'PCA', 'LDA', 'MDS', 'ProbPCA', 'Isomap', 'KernelPCA', 'GDA', 'tSNE', 'ManifoldChart'};
dimReduxMethod = 'tSNE';
YY = compute_mapping(observation, dimReduxMethod, 2, 5, 50);

%%
[~, ~, segmentMap] = models_LoadCfg(d);cfg.express = segmentMap{1};cfg.categorize = segmentMap{2};cfg.featureMin = segmentMap{3};cfg.featureMax = segmentMap{4};
getCoordinates = @(x) feval(d.categorize, observation, cfg);
mapLinIndx= getBins( observation, getCoordinates, d, edges);
xlin = -20:1:20; ylin = xlin;
[xq,yq] = meshgrid(xlin, ylin);

binU = griddata(YY(:,1),YY(:,2),double(binColorU(mapLinIndx)),xq,yq, 'natural');
binV = griddata(YY(:,1),YY(:,2),double(binColorV(mapLinIndx)),xq,yq, 'natural');
binW = griddata(YY(:,1),YY(:,2),double(binColorW(mapLinIndx)),xq,yq, 'natural');

fig(end+1) = figure(101);hold off;
imagesc(cat(3,flipud(binU),flipud(binV),flipud(binW)));
axis equal;axis tight;
title('Bin IDs');
drawnow;
