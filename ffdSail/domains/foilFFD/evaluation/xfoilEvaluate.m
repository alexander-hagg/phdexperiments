function [cD,cL] = xfoilEvaluate(coord, foamDir)

MACH = 0.5;
AOA  = 2.7;
RE = 1e6;

oldPath = pwd;
[cD, cL] = xfoilCdCl(coord', foamDir, AOA, RE, MACH,'pane oper iter 100');
eval(['cd ' oldPath])

end