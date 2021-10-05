function [map,configs,stats] = poem(observations,phenotypes,fitness,p,d)
%POEM - Phenotype Optimization through Evolution of a Manifold
% Main run script of POEM algorithm
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2019; Last revision: 12-Nov-2019

nSamples = size(observations,1);

if p.display.illu
    figure(1);hold off;axMapFit = gca;
end

for iter=1:p.numIterations
    
    % 1. Train AE
    p.map.model = trainFeatures(phenotypes,p.model);
    
    % 2. Illuminate with QD
    map = createMap(d, p.map,p.map.extraMapValues);
    features = p.categorize(observations,phenotypes,p.map,d);
    p.map.categorize = p.categorize;
    
    configs{iter} = p;
    
    [replaced, replacement] = nicheCompete(observations,fitness,map,d,p.map,features);
    map = updateMap(replaced,replacement,map,fitness,observations,features);
    
    map = illuminate(map,d.fitfun,p.map,d,p.model);
    
    % 3. Select new model members
    candidates = reshape(map.genes,[],d.dof);
    candidates(all(isnan(candidates)'),:) = [];
    observations = candidates;
    [fitness,phenotypes] = d.fitfun(observations);
    
    % 4. Statistics
    stats.fitness.mean(iter) = nanmean(map.fitness(:)); stats.fitness.median(iter) = nanmedian(map.fitness(:)); stats.fitness.std(iter) = nanstd(map.fitness(:)); stats.fitness.total(iter) = nansum(map.fitness(:));
    stats.elites.number(iter) = sum(~isnan(map.fitness(:)));
    stats.model.trainingLosses(iter,:) = p.map.model.losses;
    stats.maps{iter} = map;
    
    if p.display.illu
        viewMap(map,'fitness', d, axMapFit);caxis(axMapFit,[0 1]);drawnow;
    end
    
end

end
