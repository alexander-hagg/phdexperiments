basemesh = stlread('ffdSail/domains/mirror/ffd/mirrorBase.stl');
basemesh.vertices(:,1:2) = basemesh.vertices(:,1:2) - mean(basemesh.vertices(:,1:2));

figure(1);scatter3(basemesh.vertices(:,1),basemesh.vertices(:,2),basemesh.vertices(:,3));

stlwrite('mirrorBase.stl', basemesh.faces', basemesh.vertices');