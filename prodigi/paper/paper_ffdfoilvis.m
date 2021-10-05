mutation = 0.5*ones(1,10);
newMeshPoints =  ffdRaeY_single(mutation);
clf;
fig(1) = figure(1);hold off;
plot(newMeshPoints(1,:),newMeshPoints(2,:),'k','LineWidth', 4);
hold on;
ax = gca;
ax.Visible = 'off'
xlabel('X');ylabel('Y');

ctlX = 0:1/6:1;
ctlX = repelem(ctlX(2:end-1),1,2);
ctlY = 0.2*repmat([1 -1],1,5);
scatter(ctlX,ctlY,64,'k','filled');

for i=1:length(ctlX)
    text(ctlX(i)+0.03,ctlY(i)+0.03,int2str(i),'FontSize',16);
end

ctlX = 0:1/6:1;
ctlX = [ctlX fliplr(ctlX) 0];
ctlY = 0.2*[repelem([1 -1],1,7) 1];
plot(ctlX,ctlY,'k-', 'LineWidth',2);
%annotation('textarrow',x,y,'String','y = x ')

%axis equal;grid on;axis tight;axis([0 1 -0.25 0.25]);

save_figures(fig, './', ['ffdfoil_rep_'], 16, [7 6]);
