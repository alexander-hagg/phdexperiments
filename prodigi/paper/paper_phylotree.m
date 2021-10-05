clear selProt1 selProt2;
selProt1 = nan(6,51);selProt2 = nan(6,51);
for SEL=1:2
    for i=1:3
        % Get selected prototype
        prototypes = prodigi_output{SEL}{i}.prototypes;
        selProtID = prodigi_output{SEL}{i}.conceptSelection.id;
        
        
        if i==1
            % Get prototypes with largest distance from selected
            distances = triu(pdist2(prototypes(selProtID,:),prototypes));
            [sortDist,sortDistIDs] = sort(distances(:),'descend');   
            sortDistIDs = sortDistIDs';
        else
            % Get largest prototypes
            sizes = prodigi_output{SEL}{i}.concepts.sizes;
            sizes = sizes(2:end);
            [sortDist,sortDistIDs] = sort(sizes,'descend');
        end
        
        IDs = [selProtID sortDistIDs];
        IDs = unique(IDs,'stable');
        IDs = IDs(1:3);
        if SEL==1
            selProt1((i-1)*3+1:(i)*3,:) = prototypes(IDs,:);
        elseif SEL==2
            selProt2((i-1)*3+1:(i)*3,:) = prototypes(IDs,:);
        end
    end
end

%%
clear fig;
set(0,'DefaultFigureWindowStyle','default')

clear fig
zAngle = 107;yAngle = 0;fontsize = 16;
fig = viewMirrors(d, [selProt1;selProt2], yAngle, fontsize, zAngle, 0, false);

for i=1:length(fig)
    figure(i);
    title(int2str(i));
    ax = gca;
    ax.Title.Visible = 'on';
    saveas(gcf,['mirror_tree_' int2str(yAngle) '_' int2str(i) '.eps'],'epsc')
end
