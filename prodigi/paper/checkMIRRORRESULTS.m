set(0,'DefaultFigureWindowStyle','default')

for i=1:3
    XP = prodigi_output{1}{i};
    size(XP.prototypes,1)
    XP = prodigi_output{2}{i};
    size(XP.prototypes,1)
    numP = size(XP.prototypes,1);
    fig((i-1)*2+1) = figure((i-1)*2+1);hold off;
    scatter(XP.latent(:,1),XP.latent(:,2),'filled');
    hold on;
    [~,id] = ismember(XP.prototypes(1:numP,:),XP.optima,'rows');
    scatter(XP.latent(id,1),XP.latent(id,2),'r','filled');
    %
    
    fig((i-1)*2+2) = figure((i-1)*2+2); 
    [~,~,cLabel] = viewMap(XP.sail.predMap(end).dragForce,d);caxis([1 3]);
    colormap(parula(6));
    cLabel.Label.String = 'Drag Force [N]';
end

%save_figures(fig, './', ['mirror_maps_'], 12, [6 6]);


%% Prototypes
%close all;
set(0,'DefaultFigureWindowStyle','docked')
for i=1:3
    XP = prodigi_output{2}{i};
    XP.conceptSelection.id
    viewMirrors(d, XP.prototypes(XP.conceptSelection.id,:), [0], 16, 107,i,false);
    %viewMirrors(d, XP.prototypes(1:5,:), [90], 16, 107,2,false);    
end

%%
set(0,'DefaultFigureWindowStyle','docked')

for i=1:16
    optimum = XP.sail.predMap(end).genes(i,16,:);
    %optimum = XP.sail.predMap(end).genes(16,i,:);
    viewMirrors(d, optimum, [0], 16, 107, 0, false);
    viewMirrors(d, optimum, [90], 16, 107, 1, false);drawnow;
end

