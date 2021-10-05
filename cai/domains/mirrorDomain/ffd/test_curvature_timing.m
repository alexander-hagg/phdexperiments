%tic;Cmean=patchcurvature(FV);toc
tic;
[gm samc] = mcurvature_vec(FV.vertices(:,1)',FV.vertices(:,2)',FV.vertices(:,3)')
toc;