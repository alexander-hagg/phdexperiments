

function [y] = welchetal92(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% WELCH ET AL. (1992) FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, ..., x20]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term1 = 5.*xx(:,12) ./ (1+xx(:,1));
term2 = 5 * (xx(:,1)-xx(:,20)).^2;
term3 = xx(:,5) + 40*xx(:,19).^3 - 5*xx(:,19);
term4 = 0.05*xx(:,2) + 0.08*xx(:,3) - 0.03*xx(:,6);
term5 = 0.03*xx(:,7) - 0.09*xx(:,9) - 0.01*xx(:,10);
term6 = -0.07*xx(:,11) + 0.25*xx(:,13).^2 - 0.04*xx(:,14);
term7 = 0.06*xx(:,15) - 0.01*xx(:,17) - 0.03*xx(:,18);

y = term1 + term2 + term3 + term4 + term5 + term6 + term7;

end

