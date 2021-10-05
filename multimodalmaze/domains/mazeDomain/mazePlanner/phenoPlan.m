function phenotype = phenoPlan(genome,d)
%GENOMEPLAN Summary of this function goes here
%   Detailed explanation goes here

planLength  = length(genome)/2;
start       = d.startPosition';

phenotype = [start start+[genome(1:planLength);genome(planLength+1:end)]];



% planLength  = length(genome)/2;
% start       = d.center';
% theta       = (d.theta/180) * pi;
% 
% heading = theta;
% phenotype(:,1) = [start;heading];
% for p=1:planLength-1
%     a = pi * 2 * genome(p)/range(d.ranges);
%     l = (1 + genome(p+planLength));
%     
%     heading = heading + a;
%     heading = mod(heading,pi);
%     dx = sin(heading) * l;
%     dy = cos(heading) * l;
%     
%     phenotype(:,p+1) = phenotype(:,p) + [dx;dy;heading];
% end



end

