%
currentIteration = 0;
set(0,'DefaultFigureWindowStyle','docked')
rerunPredict = true;
    
while(true)
    %%
    system(['scp ahagg2s@wr0.wr.inf.h-brs.de:/home/ahagg2s/archetype/allMapsSPHEN.mat .']);
    
    load('allMapsSPHEN.mat');
    nMaps = length(allMaps);
    
    for i=2:nMaps
        diffMap{i} = allMaps{i};
        diffMap{i}.fitness = (allMaps{i}.fitness-allMaps{i-1}.fitness).*~isnan(allMaps{i-1}.fitness);
    end
    
    disp(['Iteration: ' int2str(nMaps)]);
    if nMaps <= currentIteration
        disp('Waiting for update')
        pause(120);
    else
        
        currentIteration = nMaps;
        disp(['Iteration ' int2str(currentIteration) ' received on: ' char(datetime('now'))])
        if rerunPredict
            disp('Rerunning prediction map creation one last time including the actual locations of samples in feature space, not predicted ones');
            pp = sphenP; %pp.display.illu = true; %pp.display.illuMod = 25;
            pp.nGens = 4096;
            [predMap, ~] = createPredictionMap([fitModel,featModels{:}],pp,d);
        end
        
        
        %
        fig(1) = figure(1);hold off;
        scatter(1:64,trueFitness(1:64),32,'k','filled');
        hold on;
        plot([1,length(trueFitness)],[mean(trueFitness(1:64)) mean(trueFitness(1:64))],'k--','LineWidth',1);
        cmap = hsv(256);
        for n=65:16:length(trueFitness)
            rndClr = randi(256);
            scatter(n:n+15,trueFitness(n:n+15),32,repmat(cmap(rndClr,:),16,1),'filled');
            plot([n,n+15],[median(trueFitness(n:n+15)) median(trueFitness(n:n+15))],'Color','k','LineWidth',2);
        end
        ax = gca;
        %ax.YAxis.Limits = [0.45 0.8];
        ax.YAxis.Limits = [0 1];
        grid on;
        ylabel('Acquisition Fitness');
        xlabel('Sample ID');
        legend('init','median fit init','acq','mean fit acq','Location','NorthWest');
        title(['Iteration: ' int2str(nMaps)]);
        drawnow;
        
        
        %%
        fig(2) = figure(2);hold off;
        %scatter(1:length(trueFitness),trueFitness,32,'k','filled');
        %hold on;
        predMapFitVals = []; filled = [];
        mapSize = (numel(allMaps{i}.fitness));
        for i=1:nMaps
            predMapFitVals(i) = nanmean(allMaps{i}.fitness(:));
            predMapFitValsMedian(i) = nanmedian(allMaps{i}.fitness(:));
            filled(i) = sum(~isnan(allMaps{i}.fitness(:))) ./ mapSize;
            maxfit(i) = nanmax(allMaps{i}.fitness(:));            
        end
        plot(80:16:64+currentIteration*16,predMapFitVals,'LineWidth',2);
        hold on;
        plot(80:16:64+currentIteration*16,filled,'LineWidth',2);
        plot(80:16:64+currentIteration*16,maxfit,'LineWidth',2);
        
        if rerunPredict
            scatter(64+currentIteration*16,nanmedian(predMap.fitness(:)),32,'b','filled')
            scatter(64+currentIteration*16,sum(~isnan(predMap.fitness(:))) ./ mapSize,32,'r','filled')
            scatter(64+currentIteration*16,nanmax(predMap.fitness(:)),32,'y','filled')
        end
        
        ylabel('Mean Prediction Fitness');
        ax = gca;
        ax.YAxis.Limits = [0 1];grid on;
        legend('median prediction fitness', 'filled % prediction','max fit','medpred 4096','filled 4096','max fit 4096','Location','NorthWest');
        drawnow;
        
        %%
        fig(3) = figure(3);
        improved = []; degraded = []; totalfit = []; maxfit = [];
        for i=1:nMaps
            if i > 1
                allFitness = diffMap{i}.fitness(:);
                improved(i) = nanmean(allFitness(allFitness>0));
                degraded(i) = abs(nanmean(allFitness(allFitness<0)));
            end
            totalfit(i) = nansum(allMaps{i}.fitness(:));
        end
        yyaxis left;hold off;
        plot(80:16:64+currentIteration*16,improved,'LineWidth',2);hold on;
        plot(80:16:64+currentIteration*16,degraded,'LineWidth',2);
        yyaxis right;hold off;
        plot(80:16:64+currentIteration*16,totalfit,'LineWidth',2);hold on;
        if rerunPredict
            scatter(64+currentIteration*16,nansum(predMap.fitness(:)),32,'r','filled')
        end
        ylabel('Total Fitness');
        ax = gca;
        ax.YAxis(1).Limits = [0 1.2*max([improved,degraded])];grid on;
        ax.YAxis(2).Limits = [0 1.2*max([totalfit maxfit])];grid on;
        legend('Avg. Improvement','Avg. Degradation','Total Fitness','Total Fit 4096','Location','NorthWest');
        drawnow;
        
        %%
        fig(4) = figure(4);
        step = ceil((nMaps-1)/7);
        
        t=1;
        for i=[1:step:(nMaps-1),nMaps]
            subplot(3,3,t);
            t=t+1;
            %%
            hold off;
            [~,~,cb] = viewMap(allACQMaps{i},d);
            hold on;
            cb.Label.String = 'Acquisition';
            caxis([0 1]);
            colormap(parula(8));
            title(['A it. ' int2str(i)]);
            
            clear tf;
            tf(:,1) = featModels{1}.trainOutput([1:64+(i-1)*16],:);
            tf(:,2) = featModels{2}.trainOutput([1:64+(i-1)*16],:);
            %scatter(32*tf(:,1)+0.5,32*tf(:,2)+0.5,8,'k','filled')
            scatter(32*tf(:,1)-0.5,32*tf(:,2)+0.5,8,'k','filled')
            
        end
        
        %%
        fig(5) = figure(5);
        step = ceil((nMaps-1)/7);       
        t=1;
        for i=[1:step:(nMaps-1),nMaps]
            subplot(3,3,t);
            t=t+1;
            %%
            [~,~,cb] = viewMap(allMaps{i},d);
            cb.Label.String = 'Pred. Fitness';
            colaxisrange = [0 1];
            caxis(colaxisrange);
            colormap(gca,parula(16));
            title(['P it. ' int2str(i)]);
        end
        subplot(3,3,9);
        %%
        if rerunPredict
            [~,~,cb] = viewMap(predMap,d);
            cb.Label.String = 'Pred. Fitness';
            colaxisrange = [0 1];
            caxis(colaxisrange);
            colormap(gca,parula(16));
            title(['4096 gen final']);
        end
        %%
        % Acquired Samples
        %featMultiplier = 50;
        figHandle = figure(7); hold off;
        fig(6) = figHandle;
        genes = featModels{1}.trainInput;
        features = [featModels{1}.trainOutput,featModels{2}.trainOutput];
        %genes = reshape(allACQMaps{end}.genes,[],d.dof);
        tRes = size(allACQMaps{end}.fitness,1);
        %x = 1:tRes; y = 1:tRes;
        %[X,Y] = ndgrid(x,y); coordinates = [X(:),Y(:)];
        coordinates = 32*features;
        coordinates(:,2) = 33-coordinates(:,2);
        
        nans = all(isnan(genes)');
        genes(nans,:) = []; coordinates(nans,:) = [];
        figHandle = showPhenotype(genes,d,figHandle,coordinates);
        axisLimits = [-0.5 tRes+0.5 -0.5 tRes+0.5];
        axis(axisLimits);
        axis equal;
        set(gca, 'Visible', 'off')
        title('Acquired Samples');
        %%
        %set(0,'DefaultFigureWindowStyle','default')
        if rerunPredict
            genes = reshape(predMap.genes,[],d.dof);
            tRes = size(allMaps{end}.fitness,1);
            
            x = 1:tRes; y = 1:tRes;
            [X,Y] = ndgrid(x,y); coordinates = [X(:),Y(:)];
            coordinates(:,2) = 33-coordinates(:,2);
            nans = all(isnan(genes)');
            genes(nans,:) = []; coordinates(nans,:) = [];
            
            figHandle = figure(8); hold off;
            fig(7) = figHandle;
            figHandle = showPhenotype(genes,d,figHandle,coordinates);
            axisLimits = [-0.5 tRes+0.5 -0.5 tRes+0.5];
            axis(axisLimits);
            axis equal;
            set(gca, 'Visible', 'off')
            
            figHandle = figure(9); hold off;
            fig(8) = figHandle;
            figHandle = showPhenotypeBMP(genes,d,figHandle,coordinates);
        end
        
        %%
        
        figHandle = figure(10); hold off;
        fig(9) = figHandle;
        genes = reshape(allMaps{end}.genes,[],d.dof);
        %genes = reshape(predMap.genes,[],d.dof);
        tRes = size(allMaps{end}.fitness,1);
        
        x = 1:tRes; y = 1:tRes;
        [X,Y] = ndgrid(x,y); coordinates = [X(:),Y(:)];
        coordinates(:,2) = 33-coordinates(:,2);
        nans = all(isnan(genes)');
        genes(nans,:) = []; coordinates(nans,:) = [];
        figHandle = showPhenotype(genes,d,figHandle,coordinates);
        axisLimits = [-0.5 tRes+0.5 -0.5 tRes+0.5];
        axis(axisLimits);
        axis equal;
        set(gca, 'Visible', 'off')
        
        
    end
    
end

save_figures(fig, '.', ['lettuceSPHEN_vis'], 12, [5 4]);

%% Create GIFs
set(0,'DefaultFigureWindowStyle','default')
nMaps = length(allMaps);
featMultiplier = 50;
filename = 'lettuceShapeEvolution.gif';
close all; clear figHandle;
for i=[1,nMaps]%1:nMaps
    figHandle(i) = figure(i);
    figHandle(i).Renderer = 'painters';
    figHandle(i).Position = [0 0 800 800];
    hold off;delete(findall(gcf,'type','annotation'));
    
    genes = reshape(allMaps{i}.genes,[],d.dof);
    fitness = reshape(allMaps{i}.fitness,[],1);
    tRes = size(allMaps{i}.fitness,1);
    x = 1:tRes; y = 1:tRes;
    [X,Y] = ndgrid(x,y); features = [X(:),Y(:)];
    features(:,2) = 33-features(:,2);
    nans = all(isnan(genes)');
    genes(nans,:) = []; features(nans,:) = []; fitness(nans) = [];
    clrs = fitness.*[0 2 0] + (1-fitness).*[1 0 0]; clrs(clrs>1) = 1;
    figHandle(i) = showPhenotype(genes,d,figHandle(i),features,clrs);
    set(gca, 'Visible', 'off')
    axisLimits = [-0.5 tRes+0.5 -0.5 tRes+0.5];
    axis(axisLimits);
    %axis tight manual % this ensures that getframe() returns a consistent size
    %axis manual
    %axis equal;
    dim = [.8 .65 .3 .3];
    str = ['Iteration ' int2str(i) '/' int2str(64)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
end
drawnow;

%%
for i=1:nMaps
    %% Capture the plot as an image
    frame = getframe(figHandle(i));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',0,'DelayTime',0.3);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.3);
    end
end



