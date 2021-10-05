function viewMirrorsMap(predMap)
%%VIEWMIRRORSMAP shows all mirror meshes in a map
%
%

set(0,'DefaultFigureWindowStyle','default')

for i=1:16
    for j=1:16
        mutation = predMap.genes(i,j,:);
        pred_fitness = predMap.fitness(i,j,:);
        [FV, ~, ffdP] = mirror_ffd_Express(mutation, 'mirrorBase.stl');
        
        fig(1) = figure(1);
        delete(findall(gcf,'type','annotation'))
        subplot(1,2,1);
        hold off;
        hm = mirrorVisPaper(FV,ffdP, false, false, false);
        hm.EdgeColor = [0 0 0];
        hm.FaceColor = [1 1 1];
        hm.FaceAlpha = 1;
        xlabel('x');ylabel('y');zlabel('z');
        view(90,90);grid on;axis equal;
        axis([750 1050 -1150 -700 500 800]);
        ax = gca;
        ax.Visible = 'off';
        
        subplot(1,2,2);
        hold off;
        hm = mirrorVisPaper(FV,ffdP);
        hm.EdgeColor = [0 0 0];
        hm.FaceColor = [1 1 1];
        hm.FaceAlpha = 1;
        str = ['Predicted Fitness: ' num2str(pred_fitness)];
        dim = [.4 .5 .3 .3];
        a = annotation('textbox',dim,'String',str,'FitBoxToText','on');
        xlabel('x');ylabel('y');zlabel('z');
        view(90,0);grid on;axis equal;
        axis([750 1050 -1150 -700 500 800]);
        ax = gca;
        ax.Visible = 'off';
        drawnow;
    end
end

end