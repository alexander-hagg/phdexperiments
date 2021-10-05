d = mirror_Domain; % for point IDs
%% Animated gif of random mirror shapes

h = figure(1);
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for n = 1:0.5:20
    FV = mirror_ffd_Express(rand(1,51), d.FfdP);
    V = FV.vertices';F = FV.faces;
    subplot(2,1,1);
    title('Random car mirrors');
    hm = plotmesh(V,F);
    hm.FaceColor = [1 1 1];
    view(100,0);axis([-150 120 -150 130 600 800]);
    
    subplot(2,1,2);
    hm = plotmesh(V,F);
    hm.FaceColor = [1 1 1];
    view(100,90);axis([-150 120 -150 130 600 800]);
    
    drawnow
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end


%% Used for construction of features
FV = mirror_ffd_Express(0.5*ones(1,51), d.FfdP);
V = FV.vertices';F = FV.faces;


fig(1) = figure(1);
hm = plotmesh(V,F);
hm.FaceColor = [1 1 1];
view(107,90);axis([-150 120 -150 130 600 800]);

fig(2) = figure(2);
scatter3(V(:,1),V(:,2),V(:,3));
view(107,90);axis([-150 120 -150 130 600 800]);

fig(3) = figure(3);
scatter3(border1(:,1),border1(:,2),border1(:,3));
view(107,0);axis([-150 120 -150 130 600 800]);

fig(4) = figure(4);
scatter3(border2(:,1),border2(:,2),border2(:,3));
view(107,0);axis([-150 120 -150 130 600 800]);

[tf, baseSubMeshIds]=ismember(submesh,V,'rows');

%% Visualize features
[FV,~,ffdP] = mirror_ffd_Express(0.5*ones(1,51), d.FfdP);
[FVRED,~,ffdPRED] = mirror_ffd_Express(0.5*ones(1,51), 'mirrorBaseNetfabbPaperVisualization.stl');
FVRED.vertices = FVRED.vertices';
%V = FV.vertices';F = FV.faces;
%[mF,mV] = reducepatch(F,V,0.3);
%FVRED.vertices = mV; FVRED.faces = mF;
%
clear fig;
fig(1) = figure(1);hold off;
hm = mirrorVisPaper(FV,ffdP,true,false,false,FVRED);
title(['Original mesh']);
feature(1) = getTotalCurvature(FV.vertices, d);
plot3(V(d.features.curvature.ids,1),V(d.features.curvature.ids,2),V(d.features.curvature.ids,3),'r-','LineWidth',4);
[feature(2), maxVertex, minVertex, maxVertexID, minVertexID] = getRelativeLength(V', d);

%
plot3([-30 66],[-164 -112],[680 680],'b-','LineWidth',4);
hold on;
%plot3(V([maxVertexID,maxVertexID],1),V([minVertexID,maxVertexID],2),[680 680],'b:','LineWidth',4);
%plot3(V([minVertexID,minVertexID],1),V([minVertexID,maxVertexID],2),[680 680],'b:','LineWidth',4);

ax = gca;ax.Visible = 'off';axis tight;
view(90,40);
axis([-120 150 -180 180 610 770]);

save_figures(fig, './', ['mirror_features_'], 12, [6 6]);

%%
fig(2) = figure(2);hold off;
hm = mirrorVisPaper(FV,ffdP,false,false,false);
% Show feature 1
feature(1) = getTotalCurvature(FV.vertices, d);
plot3(V(d.features.curvature.ids,1),V(d.features.curvature.ids,2),V(d.features.curvature.ids,3),'r-','LineWidth',4);
title(['Total Curvature']);

ax = gca;ax.Visible = 'off';
view(107,0);
axis([-150 120 -180 120 610 770]);

%
fig(3) = figure(3);hold off;
hm = mirrorVisPaper(FV,ffdP,false,false,false);
% Show feature 2
[feature(2), maxVertex, minVertex, maxVertexID, minVertexID] = getRelativeLength(V', d);
plot3(V([minVertexID,maxVertexID],1),V([maxVertexID,maxVertexID],2),[760 760],'r-','LineWidth',4);
hold on;
plot3(V([maxVertexID,maxVertexID],1),V([minVertexID,maxVertexID],2),[760 760],'r:','LineWidth',4);
plot3(V([minVertexID,minVertexID],1),V([minVertexID,maxVertexID],2),[760 760],'r:','LineWidth',4);
title(['Length: ' num2str(feature(2))]);
view(90,90);
axis([-150 120 -180 120 610 770]);
ax = gca;ax.Visible = 'off';

save_figures(fig, './', ['mirror_features_'], 12, [10 10]);

%feature(3) = getMirrorSurface(FV.vertices, d);

%% 

%% Visualize curvature
FV = mirror_ffd_Express(0.5*ones(1,51), d.FfdP);
V = FV.vertices';F = FV.faces;

fig(6) = figure(6);
subplot(2,1,1);hold off;
hm = plotmesh(V,F);hold on;hm.FaceColor = [1 1 1];
% Show feature 1
feature(1) = getTotalCurvature(FV.vertices, d);
plot3(V(d.features.curvature.ids,1),V(d.features.curvature.ids,2),V(d.features.curvature.ids,3),'r-','LineWidth',4);
title(['Total Curvature: ' num2str(feature(1))]);
view(107,0);axis([-150 120 -150 130 600 800]);%axis tight;

mutation = 0.9*[-1 repmat([1 -1],1,25)];
mutation = rand(1,51);
FV = mirror_ffd_Express(mutation, d.FfdP);
V = FV.vertices';F = FV.faces;

subplot(2,1,2);hold off;
hm = plotmesh(V,F);hold on;hm.FaceColor = [1 1 1];
% Show feature 1
feature(1) = getTotalCurvature(FV.vertices, d);
plot3(V(d.features.curvature.ids,1),V(d.features.curvature.ids,2),V(d.features.curvature.ids,3),'r-','LineWidth',4);
title(['Total Curvature: ' num2str(feature(1))]);
view(107,0);axis([-150 120 -150 130 600 800]);%axis tight;



%%
h(1) = patch('Vertices',FV.vertices','Faces',FV.faces, 'FaceAlpha',0.9,'FaceColor', [0.7 0.7 0.7], 'LineStyle', 'none');
hold on;
% Plot curvature lines
for i=1:length(d.features.curvature.ids)
    pp = FV.vertices(:,d.features.curvature.ids{i});
    h(2) = scatter3(pp(1,:),pp(2,:),pp(3,:),32,'r','filled');
end
% Plot relative length
pp = FV.vertices(:,[minVertexID, maxVertexID]);
h(3) = scatter3(pp(1,:),pp(2,:),pp(3,:),128,'g','filled');

% Plot area
pp = FV.vertices(:,d.features.mirror.ids);
pp = [pp pp(:,1)];
h(4) = plot3(pp(1,:),pp(2,:),pp(3,:),'k','LineWidth',4);

legend([h(1) h(2) h(3) h(4)],'Mirror', 'Curvature', 'Relative Length', 'Mirror Area');
title(['Features: ' num2str(feature(1),3) ' - ' num2str(feature(2),3) ' - ' num2str(feature(3),3)]);

lightangle(-65, 20);
lighting phong; light; material metal;
view(-50,30);
grid on;
xlabel('x');ylabel('y');zlabel('z');


%save_figures(fig, './', ['feature_visualization_'], 12, [16 10]);
