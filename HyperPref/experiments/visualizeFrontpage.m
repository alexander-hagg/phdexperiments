addpath(genpath(pwd));
load('catmullRom_step4');

%% Show hypervolumes and selected shapes
for i=1:4
    [genes,fitness,features,bins] = extractMap(maps{i});
    
    
    accepted = fitness > 0.1;
    axHandle = showPhenotype(genes(accepted,:), d, 1, [], bins(accepted,:)); 
    figs(i) = axHandle.Parent;
    if i < 4
        placement = [d.selection.selected{i}(:,1),d.selection.selected{i}(:,2)];
        selbin = [];
        for iDim = 1:size(placement,2)
            selbin(:,iDim) = discretize(placement(:,iDim),maps{i}.edges{iDim});
        end
        %placement = placement * p.featureResolution(1);
        selbin(:,2) = -selbin(:,2);        
        scatter(selbin(:,1),selbin(:,2),32,[1 0 0],'filled');
    end
    
    title(['Iteration ' int2str(i)]);
end

%% Save figures

save_figures(figs, '.', 'rand2human', 16, [5 5])