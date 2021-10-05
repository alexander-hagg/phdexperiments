selPrototypes = 3;
set(0,'DefaultFigureWindowStyle','default')

fontsize = 16;
close all;clf;
for selIter=1:3
    clear fig;
    sizes = output{selIter}.concepts.sizes(2:end);
    [~,sizeSortID] = sort(sizes,'descend');sizeSortID = sizeSortID(1:selPrototypes);
    fig = viewMirrorPrototypes(output{selIter}.prototypes(sizeSortID,:),sizes(sizeSortID),[0 90],fontsize);
    save_figures(fig, './', ['mirror_ITER_' int2str(selIter) '_prototypes_'], fontsize, [5 4]);
end



