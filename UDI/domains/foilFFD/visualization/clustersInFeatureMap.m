%% Clusters in feature map
%close all;
clear fig; set(0,'DefaultFigureWindowStyle','docked');

for it=2:length(output)
    for MAINITER = 1:2
        clusterDistMap = 99999*ones(25,25);
        clusterIDMap = nan(25,25);
        selectedPrototype = output{it}.data{1}.conceptSelection.id;
        centersLatent = output{it}.data{MAINITER}.prototypesLatent;
        centers = output{it}.data{MAINITER}.prototypes;
        
        if MAINITER==1
            numClassesONE = size(centers,1);
        end
        colors = colormap(parula(numClassesONE));
        colors(selectedPrototype,:) = [1 0 0];

        for protoID=1:size(centers,1)
            centerLatent = centersLatent(protoID,:);
            center = centers(protoID,:);
            
            XXXPAR = reshape(output{it}.data{MAINITER}.sail.predMap(end).genes,625,10);
            invalidshapes = any(isnan(XXXPAR),2);
            
            XXX = output{it}.data{MAINITER}.latent;
            %dists = pdist2(center,XXX);
            dists = pdist2(centerLatent,XXX);
            distMap = nan(25,25);
            distMap(~invalidshapes) = dists;
            
            closerCells = distMap<clusterDistMap;
            clusterDistMap(closerCells) = distMap(closerCells);
            clusterIDMap(closerCells) = protoID;
        end
        
        
        fig(it-1) = figure(it-1);
        subplot(2,1,MAINITER);
        [figHandle, imageHandle, cHandle] = viewMap(clusterIDMap,d);
        cHandle.Label.String = 'Design Class Membership';
        %caxis([0 20]);
        if MAINITER==1
            colormap(gca,colors);
        else
            %colormap(gca,colors(selectedPrototype*20:(selectedPrototype+1)*20,:));
        end
        title(['Design Classes in iteration ' int2str(MAINITER)]);
    end
end
%%
for i=1:length(fig)
    figure(i)
    name = ['sail_vs_produqd_DIST_'];
    %save2pdf([name int2str(i) '.pdf'],fig(i),600,[7 6]);
    %print([name int2str(i)],'-depsc');
    saveas(gcf,['ffd_distmap_' int2str(i) '.eps'],'epsc')
end

%% Aggregated map compared to SAIL and single maps

fitMAPSAIL = output{1}.data{1}.sail.predMap(end).fitness_true;
fig(1) = figure(1);
subplot(1,3,1);
[figHandle, imageHandle, cHandle] = viewMap(fitMAPSAIL,d);
caxis([-5 -3]);
title('SAIL true fitness');
    
fitMAP = fitMAPSAIL;
fitMAP(isnan(fitMAPSAIL)) = 999;

branchMAP = zeros(25,25);
for MAINITER = 2:2
    for it=2:length(output)
        fits = output{it}.data{MAINITER}.sail.predMap(end).fitness_true;
        betterFits = (fits < fitMAP);
        fitMAP(betterFits) = fits(betterFits);
        branchMAP(betterFits) = it;
    end
    fitMAP(fitMAP==999) = nan;
    
    fig(MAINITER) = figure(MAINITER);
    subplot(1,3,1);
    [figHandle, imageHandle, cHandle] = viewMap(fitMAP,d);
    caxis([-5 -3]);
    title('PRODUQD aggregated true fitness');
    
    subplot(1,3,2);
    colors = parula(25); colors = [[0 1 0];colors];
    percImprov = -100*(fitMAPSAIL-fitMAP)./fitMAPSAIL;
    percImprov(isnan(fitMAPSAIL)) = -1;
    percImprov(isnan(fitMAP)) = nan;
    [figHandle, imageHandle, cHandle] = viewMap(percImprov,d);
    cHandle.Label.String = '% improvement';
    caxis([-1 25]);
    colormap(gca,colors);
    title('Improvement over SAIL');
    
    subplot(1,3,3);
    [figHandle, imageHandle, cHandle] = viewMap(branchMAP,d);
    cHandle.Label.String = 'Prototype ID';
    title('Improvement from which prototype?');
    caxis([0 11]);
    %colormap(gca,[[1 0 0];parula(11)]);
    colormap(gca,[[1 1 1];hsv(11)]);
    
end


