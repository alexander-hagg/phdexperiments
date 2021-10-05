function p = visFoil(foil)
%VISFOIL Summary of this function goes here
%   Detailed explanation goes here
newMeshPoints =  ffdRaeY_single(foil');
p = plot(newMeshPoints(1,:),newMeshPoints(2,:),'LineWidth',2);
axis equal;grid on; axis([0 1 -0.2 0.2]);
end

