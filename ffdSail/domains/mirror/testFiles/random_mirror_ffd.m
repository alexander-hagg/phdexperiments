numShapes = 1;
set(0,'DefaultFigureWindowStyle','docked')

for shapeID=1:numShapes
    disp(['Generating shape ' int2str(shapeID) '/' int2str(numShapes)]);
    mutation = 0.5*ones(1,48);
    limitedValueIDsCenter = [11,29,43,14,32,45,17,35,47];
    limitedValueIDsOut = [2, 20, 37, 5, 23, 39, 8, 26, 41];
    mutation(limitedValueIDsCenter) = 0.816;
    mutation(limitedValueIDsOut) = 0.091;
    %mutation = rand(1,48);    
    
    [FV, ~, ffdP] = mirror_ffd_Express(mutation, 'mirrorBase.stl');
    
    fig(shapeID) = figure(shapeID);
    hold off;
    hm = mirrorVisPaper(FV,ffdP);
    hm.EdgeColor = [0 0 0];
    hm.FaceColor = [1 1 1];
    %hm.FaceAlpha = 0;
    title(['Params median ' num2str(median(mutation))]);
    xlabel('x');ylabel('y');zlabel('z');
    view(90,45);grid on;axis equal;
    axis([750 1050 -1150 -700 500 800]);
    ax = gca;
    ax.Visible = 'off'
end

%%
%save_figures(fig, './', ['mirror_curvature_'], 12, [5 5]);
