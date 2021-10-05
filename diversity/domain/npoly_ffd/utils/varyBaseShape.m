function phenotypes = varyBaseShape(dof,baseShape,numShapes,scaling)
%VARYBASESHAPE Scales a polygon according to scaling vector. Each vector
% element is turned into a scaled version of the polygon base shape
% 
% Inputs:
%   - scaling: vector containing scaling factors


basePheno = getPhenotype(zeros(1,dof),baseShape);
minScale = scaling(1); maxScale = scaling(2); 

scalingFactors = (minScale:(maxScale-minScale)/(numShapes-1):maxScale);
% Vary scale in polar coordinates
[theta,rho] = cart2pol(basePheno{1}.Vertices(:,1),basePheno{1}.Vertices(:,2));
rho = rho.*scalingFactors;
[x,y] = pol2cart(repmat(theta,1,size(rho,2)),rho);
% Reconstruct phenotypes
zeroGenes = zeros(length(scalingFactors),dof);
phenotypes = getPhenotype(zeroGenes,baseShape);
for i=1:length(phenotypes)
    phenotypes{i}.Vertices = [x(:,i),y(:,i)];
end

end

