set(0,'DefaultFigureWindowStyle','default')
numDOF = d.dof;
mutation = 0.5*ones(1,numDOF);
limitedValueIDsCenter = [11,29,43,14,32,45,17,35,47];
limitedValueIDsOut = [2, 20, 37, 5, 23, 39, 8, 26, 41];
mutation(limitedValueIDsCenter) = 0.816;
mutation(limitedValueIDsOut) = 0.091;

[FV, ~, ffdP] = mirror_ffd_Express(mutation, 'mirrorBase.stl');

fig(1) = figure(1);
hm = mirrorVisPaper(FV,ffdP,true,true,true);
hm.EdgeColor = [0 0 0];hm.FaceColor = [1 1 1];
title(['Params median ' num2str(median(mutation))]);
xlabel('x');ylabel('y');zlabel('z');
view(90,45);grid on;axis equal;
ax = gca;ax.Visible = 'off';
%%
fig(2) = figure(2);
hm = mirrorVisPaper(FV,ffdP,false,false,false,true);
hm.EdgeColor = [0 0 0];hm.FaceColor = [1 1 1];
title(['Params median ' num2str(median(mutation))]);
xlabel('x');ylabel('y');zlabel('z');
view(-90,0);grid on;axis equal;
ax = gca;ax.Visible = 'off';

%%
fig(3) = figure(3);
hm = mirrorVisPaper(FV,ffdP,false,false,false,false,true);
hm.EdgeColor = [0 0 0];hm.FaceColor = [1 1 1];
title(['Params median ' num2str(median(mutation))]);
xlabel('x');ylabel('y');zlabel('z');
view(-90,90);grid on;axis equal;
ax = gca;ax.Visible = 'off';



%%

save_figures(fig, './', ['mirror_features_'], 12, [10 10]);
