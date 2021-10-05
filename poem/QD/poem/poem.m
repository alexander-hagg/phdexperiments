function [map,manifold,stats] = poem(observations,fitness,phenotypes,p,d,manifold)
%POEM - Phenotype Optimization through Evolution of a Manifold
% Main run script of POEM algorithm
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2019; Last revision: 12-Nov-2019

nSamples = size(observations,1);

if p.display.illu
    figure(1);hold off;axMapFitPrior = gca;
    figure(2);hold off;axMapRecErrPrior = gca;
    figure(3);hold off;axMapFit = gca;
    figure(4);hold off;axMapUncertainty= gca;
end

%iter = 1;
%while nSamples <= p.numTotalSamples

for iter=1:p.numIterations
    
    % 1. Train manifold
    manifold = manifold.train(phenotypes,manifold);
        
    % 2. Illuminate with QD
    map = createMap(d, p.map,p.map.extraMapValues);
    latent = manifold.predict(phenotypes,manifold);

    [replaced, replacement] = nicheCompete(observations,fitness,map,d,p.map,latent);
    uncertainty = manifold.uncertainty(observations,latent,manifold,d.getPhenotype);
    map = updateMap(replaced,replacement,map,fitness',observations,latent,p.map.extraMapValues,{uncertainty});
    
    if p.display.illu
        viewMap(map,'fitness',axMapFitPrior);caxis(axMapFitPrior,[0 1]);drawnow;
        viewMap(map,'reconstructionError',axMapRecErrPrior);caxis(axMapUncertainty,[0 1000]);drawnow;
    end
    
    map = illuminate(map,d.fitfun,p.map,d,manifold);
    
    % 3. Select new manifold members
    candidates = reshape(map.genes,[],d.dof);
    candidates(all(isnan(candidates)'),:) = [];
    if strcmp(p.selectionMethod,'all')
        observations = [];fitness=[];phenotypes=[];
        newobservations = selectRandom(candidates,p);
    else
        numNewSolutions = ceil((p.selectPerc/100)*size(observations,1));
        observations(1:numNewSolutions,:) = []; fitness(1:numNewSolutions) = []; phenotypes(1:numNewSolutions,:) = [];    
        if strcmp(p.selectionMethod,'maxuncertainty')
            uncertainty = reshape(map.uncertainty,[],1);
            uncertainty(isnan(uncertainty)) = [];
            newobservations = selectUncertainty(candidates,uncertainty,p,'max');
        elseif strcmp(p.selectionMethod,'minuncertainty')
            uncertainty = reshape(map.uncertainty,[],1);
            uncertainty(isnan(uncertainty)) = [];
            newobservations = selectUncertainty(candidates,uncertainty,p,'min');
        elseif strcmp(p.selectionMethod,'random')
            newobservations = selectRandom(candidates,p);
        end
    end
    [newFitness,newPhenotypes] = d.fitfun(newobservations);
    observations = [observations;newobservations]; fitness = [fitness,newFitness]; phenotypes = [phenotypes;newPhenotypes];
    
    % 4. Statistics
    stats.uncertainty.mean(iter) = nanmean(map.uncertainty(:));
    stats.uncertainty.median(iter) = nanmedian(map.uncertainty(:));
    stats.uncertainty.std(iter) = nanstd(map.uncertainty(:));
    
    stats.fitness.mean(iter) = nanmean(map.fitness(:));
    stats.fitness.median(iter) = nanmedian(map.fitness(:));
    stats.fitness.std(iter) = nanstd(map.fitness(:));
    stats.fitness.total(iter) = nansum(map.fitness(:));
    
    stats.elites.number(iter) = sum(~isnan(map.fitness(:)));
    
    %stats.model.trainingLosses(iter,:) = manifold.model.losses;
    stats.maps{iter} = map;
    
    if p.display.illu
        viewMap(map,'fitness',axMapFit);caxis(axMapFit,[0 1]);drawnow;
        viewMap(map,'uncertainty',axMapUncertainty);caxis(axMapUncertainty,[0 1000]);drawnow;
    end
    
    iter = iter + 1;
end

end
