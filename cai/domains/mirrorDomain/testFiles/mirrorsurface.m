[FV, ~, ffdP] = mirror_ffd_Express(0.5 + 0.0*rand(1,41), 'mirrorBase.stl');
% Visualize
fig(1) = figure(1);hold off;
plot3(FV.vertices(1,:), FV.vertices(2,:), FV.vertices(3,:),'x');
hold on;

%%
%fig(1) = figure(1);hold off;
%plot3(mirrorSurface(:,1), mirrorSurface(:,2), mirrorSurface(:,3),'x');
%hold on;
%[~,id] = intersect(FV.vertices', mirrorSurface,'rows')

%%
area = getMirrorSurface(FV.vertices, d)
plot3(FV.vertices(1,d.base.mirrorIDs), FV.vertices(2,d.base.mirrorIDs), FV.vertices(3,d.base.mirrorIDs),'x');
xlabel('x');ylabel('y');zlabel('z');
axis([700 1050 -1200 -700 550 800]);
view(22,0);
view(112,0);

%%
submesh = (submesh'*ffdP.rotMat)';
fig(2) = figure(2);hold off;
plot(submesh(1,:), submesh(2,:),'x');
