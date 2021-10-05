function [validInds,phenotype] = validateChildren_plan(newChildren,d)
%validateChildren
%
%
% Author: Alexander Hagg
% Bonn-Rhein-Sieg University of Applied Sciences (HBRS)
% email: alexander.hagg@h-brs.de
% Nov 2018; Last revision: 02-Nov-2018
%
%------------- BEGIN CODE --------------
%% All children are valid
validInds = true(1, size(newChildren,1));

% Calculate validity
for i=1:size(newChildren,1)
    phenotype(i,:,:) = phenoPlan(newChildren(i,:),d);
end
    
for i=1:size(newChildren,1)
    [cx,cy,c] = improfile(d.map,phenotype(i,1,:),phenotype(i,2,:),d.phenotypeLength);    
    %d.phenotypeLength
    cx = cx(~isnan(c));
    cy = cy(~isnan(c));
    c = c(~isnan(c));
    
    validInds(i) = ~logical(sum(~c));
    
end

end

