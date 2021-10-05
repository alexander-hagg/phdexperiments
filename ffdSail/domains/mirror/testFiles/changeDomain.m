d = mirror_Domain('hpc', false, 'lowres', false, 'mirrorOnly', true, 'nCases', 4, 'runFolder', '/scratch/ahagg2s/sailCFD');

%% Check base mesh
fname = 'mirrorBase.stl';
rawStl = stlread(fname);
[STLMeshpoints,faces] = patchslim(rawStl.vertices, rawStl.faces);

% Normalize mesh points
[STLMeshpoints,normalizationFactors] = mapminmax(STLMeshpoints',-1,1);
STLMeshpoints = STLMeshpoints';

% Rotate mesh to align with bounding box
theta           = atan((0.88-0.45)/1.28);
rotMat          = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
STLMeshpoints   = STLMeshpoints * rotMat;

% Define bounding box
%maxY = 0.4;
maxY = 0.61;

% Select mesh points within bounding box
submesh         = (STLMeshpoints(:,2) + 0.37*STLMeshpoints(:,1))< maxY;
meshPoints      = STLMeshpoints(submesh,:);
nMeshPoints     = size(meshPoints,1);

figure(3);scatter3(STLMeshpoints(:,1),STLMeshpoints(:,2),STLMeshpoints(:,3));
ax = gca;
ax.YTick = [-1.5:0.1:1];
view(0,90)
figure(4);scatter3(meshPoints(:,1),meshPoints(:,2),meshPoints(:,3));
view(0,90)

%%
d = mirror_Domain('hpc', false, 'lowres', false, 'mirrorOnly', true, 'nCases', 4, 'runFolder', '/scratch/ahagg2s/sailCFD');


mutation = 1*ones(1,d.dof);
[FV, ~, ffdP] = mirror_ffd_Express(mutation, 'mirrorBase.stl');

figure(1);clf;
hm = mirrorVisPaper(FV,ffdP, true, true);
xlabel('x');ylabel('y');zlabel('z');
view(283,90)
title(['Mean mutation: ' num2str(mean(mutation)-0.5)]);
axis equal;

mutation = 0*ones(1,d.dof);
[FV, ~, ffdP] = mirror_ffd_Express(mutation, 'mirrorBase.stl');

figure(2);clf;
hm = mirrorVisPaper(FV,ffdP, true, true);
xlabel('x');ylabel('y');zlabel('z');
view(283,90)
title(['Mean mutation: ' num2str(mean(mutation)-0.5)]);
axis equal;

%%

mutation = rand(1,d.dof);
[FV, ~, ffdP] = mirror_ffd_Express(mutation, 'mirrorBase.stl');

figure(2);clf;
hm = mirrorVisPaper(FV,ffdP, true, true);
xlabel('x');ylabel('y');zlabel('z');
%view(283,90)
view(103,90)
title(['Mean mutation: ' num2str(mean(mutation)-0.5)]);
axis equal;
axis([800 1050 -1100 -700 550 800]);



