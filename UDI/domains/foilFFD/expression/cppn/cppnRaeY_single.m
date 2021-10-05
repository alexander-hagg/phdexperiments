function [ foil, flatAct ] = cppnRaeY( wVector)
%CPPNRAEY Converts a 10 dimensional weight vector to a 10 pt transform
% CPPN is structured with a 1 dimensional substrate and 10 weights
%
%       x_pos
%      /     \      [2 weights]
%   (sin)   (gaus)
%     |   X   |     [4 weights]
%   (sin)   (gaus)
%     |   X   |     [4 weights]
%  (tanh)   (tanh)
%     |       |     [0 weights]
%   y_top    y_bot
%
%
%%


substrate = linspace(0,1,5); % Static for now
w1 = nan(1,2);
w2 = nan(2,2);
w3 = nan(2,2);

% Assign Weights
weightStrength = 1;
wVector = weightStrength .* (wVector*2)-1;     % Genome between 0 and 1
w1(:) = wVector([1:2]);
w2(:) = wVector([3:6]);
w3(:) = wVector([7:10]);

% First Layer [NX1] * [1X2] = [NX2]
in = substrate;
act1  = in'*w1;
act1(:,1) = mysine(act1(:,1));
act1(:,2) = mygaus(act1(:,2));

% Second Layer [NX2] * [2X2] = [NX2]
act2 = act1*w2;
act2(:,1) = mysine(act2(:,1));
act2(:,2) = mygaus(act2(:,2));

% Third Layer [NX2] * [2X2] = [NX2]
act3 = act2*w3;
act3(:,1) = mytanh(act3(:,1));
act3(:,2) = mytanh(act3(:,2));

a1 = act1;a2=act2;a3=act3;
flatAct = act3(:)'; % Interleave vectors
foil = ffdRaeY(flatAct);

end

