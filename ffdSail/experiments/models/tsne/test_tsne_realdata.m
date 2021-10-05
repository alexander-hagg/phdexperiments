clear fig;
set(0,'DefaultFigureWindowStyle','default')
no_dims = 2; initial_dims = 2; perplexity= 100; theta = 0.5; alg = 'svd'; max_iter = 1000;

perplexities = [20 50 100 200];
for pp=1:length(perplexities)
    for rep=1:3
        % t-SNE
        reducedParSpace(pp,rep,:,:) = fast_tsne(in, no_dims, initial_dims, perplexities(pp), theta, alg, max_iter);
    end
end

%%

for pp=1:length(perplexities)
%for pp=length(perplexities)
%    for rep=1
    for rep=1:3
        disp(['Clustering perplexity: ' int2str(pp) ', rep: ' int2str(rep)]);
        distances = pdist2(squeeze(reducedParSpace(pp,rep,:,:)),squeeze(reducedParSpace(pp,rep,:,:)));
        mn = mean(distances(:));
        epsilon(pp,rep) = 0.5*mn/18.6;
        coreneighbours = 3;
        [~,cLab,cen] = dbscan(squeeze(reducedParSpace(pp,rep,:,:))', epsilon(pp,rep), coreneighbours);
        clustLabel{pp,rep} = cLab;
        centers{pp,rep} = cen;
    end
end

%%
set(0,'DefaultFigureWindowStyle','default')
        
for pp=1:length(perplexities)
%for pp=length(perplexities)
%    for rep=1
    for rep=1:3
        %% Density clustering
        disp(['------']);
        disp([int2str(max(clustLabel{pp,rep}))  ' clusters']);
        disp([int2str(sum(clustLabel{pp,rep}==0)) '/' int2str(size(reducedParSpace,3))  ' noise points']);
        
        % Visualization
        fig((pp-1)*3 + rep) = figure((pp-1)*3 + rep);hold off;
        cmap = colorcube(ceil(max(clustLabel{pp,rep}))+2);
        cmap(end,:) = []; % Remove white
        clrs = cmap(clustLabel{pp,rep}+1,:);
        clrs(clustLabel{pp,rep}==0,:) = repmat([0 0 0],sum(clustLabel{pp,rep}==0),1);
        
        for i=1:length(selected)
            numclu(i) = 100*length(unique(clustLabel{pp,rep}(pickedSamples{i})))/length(unique(clustLabel{pp,rep}));
            disp(['Method found ' num2str(numclu(i)) '% of total clusters']);
        end
        scatter(reducedParSpace(pp,rep,:,1),reducedParSpace(pp,rep,:,2),8,clrs,'filled');
        %scatter(reducedParSpace(pp,rep,pickedSamples{i},1),reducedParSpace(pp,rep,pickedSamples{i},2),8,'filled');
        axis equal;axis tight;hold on;            
        values = reducedParSpace(pp,:,:);
        axis([min(values(:)) max(values(:)) min(values(:)) max(values(:))]);
        ax = gca;
        grid on;
        %ax.XTickLabel = [];ax.YTickLabel = [];ax.ZTickLabel = [];
        ax.XTick(~mod(ax.XTick,5) & mod(ax.XTick,10)) = [];
        ax.YTick(~mod(ax.YTick,5) & mod(ax.YTick,10)) = [];
        title(['Optima found (various t-SNE and clustering runs)']);
        str = {['t-SNE perplexity: ' int2str(perplexities(pp))], ['DBSCAN eps: ' num2str(round(epsilon(pp,rep),2))], ...
            [int2str(max(clustLabel{pp,rep}))  ' clusters'], [int2str(sum(clustLabel{pp,rep}==0)) '/' int2str(size(reducedParSpace,3))  ' unassigned'], ...
            ['BN found ' num2str(round(mean(numclu(1:4)))) '% clusters'], ['GP found ' num2str(round(mean(numclu(5:8)))) '% clusters']};
        
        annotation('textbox',[.67 .6 .3 .3],'String',str,'FitBoxToText','on');
        drawnow;
        %%
    end
end
save_figures(fig, './', ['various_tsne'], 14, [7 7]);