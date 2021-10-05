set(0,'DefaultFigureWindowStyle','default');

clear normFcn;
close all;
clrs = parula(5);

for ff=1:2
    fig(ff) = figure(ff);hold off;
    
    for i=1:size(allresults,1)
        %subplot(4,1,i);
        for s=1:numMaps
            for j=1:size(allresults,2) % Replicates
                conf = (allresults{i,j}.predMap(s*50).confidence);
                CONFMAP(s,(j-1)*length(conf(:))+1:(j)*length(conf(:))) = conf(:)';
            end
        end
        %
        if ff==2
            normCONFMAP = CONFMAP;
            maxConf = max(CONFMAP');
            for ii = 1:19
                normCONFMAP(ii+1,:) = CONFMAP(ii+1,:)/maxConf(ii);
            end
            CONFMAP = normCONFMAP;
        end
        le = prctile(CONFMAP,25,2)';
        ri = prctile(CONFMAP,75,2)';
        patch([1:numMaps,numMaps:-1:1],[ri fliplr(le)], clrs(i,:),'FaceAlpha',.3, 'EdgeColor', clrs(i,:), 'LineWidth', 2);
        hold on;
        h(i) = plot(nanmedian(CONFMAP,2),'LineWidth',2,'Color',clrs(i,:));
        grid on;
        
        drawnow;
    end
    legend(h, strrep(modelNames,'_','\_'),'Location','NorthWest');
    title('Median Confidence');
    ylabel('Var');
    xlabel('Samples');
    ax = gca;
    ax.XTick = 0:10:numMaps;
    ax.XTickLabel = {string(ax.XTick*50)};
    axis([0 20 0 0.8]);    
end

save_figures(fig, './', 'variance_output_', 32, [16 8]);


