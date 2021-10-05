[clustLabel, varType] = dbscan(in(1:2000,:), 50, 0.5);
%%
figure(101);hold off;
uniques = unique(clustLabel);
for j=1:length(uniques)
    i = uniques(j);
    scatter(reducedParSpace(clustLabel==i,1),reducedParSpace(clustLabel==i,2), 'filled');
    hold on; 
end
if length(uniques) > 1
    legend(['noise'; string(uniques(2:end))]);
else
    legend('noise');
end