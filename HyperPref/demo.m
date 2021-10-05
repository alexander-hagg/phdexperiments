clear;clc;
%% Configuration
addpath(genpath(pwd))                           % Set path to all modules
DOF = 16;DOMAIN = 'catmullRom';                 % Degrees of freedom, Catmull-Rom spline domain
d = domain(DOF);                                % Domain configuration
p = defaultParamSet;                            % Base Quality Diversity (QD) configuration (MAP-Elites)
m = cfgLatentModel('data/workdir',d.resolution);% VAE configuration
pm = poemParamSet(p,m);                    % Configure POEM ("Phenotypic niching based Optimization by Evolving a Manifold")
pm.categorize = @(geno,pheno,p,d) predictFeatures(pheno,p.model);  % Anonymous function ptr to phenotypic categorization function (= VAE)
hypID = 1;

%% Initialize Experiment
% Initialize solution set using space filling Sobol sequence in genetic space
sobSequence = scramble(sobolset(d.dof,'Skip',1e3),'MatousekAffineOwen');  sobPoint = 1;
initSamples = range(d.ranges').*sobSequence(sobPoint:(sobPoint+pm.map.numInitSamples)-1,:)+d.ranges(:,1)';
[fitness,polygons] = fitfun(initSamples,d);

%% Run POEM's first iteration
[map{hypID}, config{hypID}, results{hypID}] = poem(initSamples,polygons,fitness,pm,d,2);
save([DOMAIN '_step1.mat']);
disp('HyperPref Step 1 Done');

%% Reload and extract results of first iteration and select IDs of shapes
load([DOMAIN '_step1.mat']);
[genes,fitness,features,bins] = extractMap(results{1}.maps{1});
% Visualization
colors = [zeros(size(genes,1),1) fitness zeros(size(genes,1),1)];            
showPhenotype(genes, d, 1.1, [], bins, colors); title('1st Iteration Result');

%% Selection
currentSelection = randi(numel(fitness)); % In this demo, select one random shape
d.selection.models{hypID} = results{hypID}.models; % Save user model to use as constraint model

phenotypes = d.getPhenotype(genes);
d.selection.selected{hypID} = features(currentSelection,:); 
d.selection.deselected{hypID} = features; d.selection.deselected{hypID}(currentSelection,:) = [];

hypID = hypID + 1;
newSamples = genes(currentSelection,:);
nNewPerSelected = ceil(pm.map.numInitSamples./length(currentSelection));
% Perturb selected shapes
for i=1:length(currentSelection)
    newSampleMutations = pm.mutSelection * randn(nNewPerSelected,d.dof);
    newSamples = [newSamples; genes(currentSelection(i),:) + newSampleMutations];
end
[newSamplesfitness,newSamplespolygons] = fitfun(newSamples,d); % Recalculate fitness! (User selection influences fitness values)
            
% figure(99);plot(sort(newSamplesfitness));hold on;
% showPhenotype(newSamples,d,p.featureResolution(1),[]); title('Injected Perturbations of Selection');

%% Run POEM's second iteration based on the user selection
[map{hypID}, config{hypID}, results{hypID}] = poem(newSamples,newSamplespolygons,newSamplesfitness,pm,d,2);
save([DOMAIN '_step2.mat']);
disp('HyperPref Step 2 Done');

%% Reload and extract results of third iteration and visualize
load([DOMAIN '_step2.mat']);
[genes,fitness,features,bins] = extractMap(results{2}.maps{1});
% Visualization
colors = [zeros(size(genes,1),1) fitness zeros(size(genes,1),1)];            
showPhenotype(genes, d, p.featureResolution(1), [], bins, colors); title('1st Iteration Result');

