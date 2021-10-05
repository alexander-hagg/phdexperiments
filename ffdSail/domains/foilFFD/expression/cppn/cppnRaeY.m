function [ foil, ffdGenome ] = cppnRaeY( wVector)
%CPPNRAEY Converts a 10 dimensional weight vector to a 10 pt transform
%
%   wVector [NX10] between 0 and 1
%
%
% CPPN is structured with a 1 dimensional substrate and 10 weights
%
%       x_pos
%      /     \      [2 weights]
%   (sin)   (gaus)      a1
%     |   X   |     [4 weights]
%   (sin)   (gaus)      a2
%     |   X   |     [4 weights]
%  (tanh)   (tanh)      a3
%     |       |     [0 weights]
%   y_top    y_bot   ffdGenome
%
%
%%
substrate = linspace(0,1,5)'; % Static for now
nGenomes  = size(wVector,1);
ffdGenome = nan([size(wVector,1) length(substrate)*2]);    % Will change with larger substrate

% Assign Weights
weightStrength = 1;
wVector = weightStrength .* (wVector*2)-1;     % Genome between 0 and 1, Weights between -1 and 1
w1 = reshape( wVector(:,1:2 )', [1 2 nGenomes]);
w2 = reshape( wVector(:,3:6 )', [2 2 nGenomes]);
w3 = reshape( wVector(:,7:10)', [2 2 nGenomes]);

in = substrate;
for i=1:nGenomes % TODO: must be a way to vectorize this
    a1 = in .* w1(:,:,i);
    a1 = sinGaus(a1');
    a2 = a1 *  w2(:,:,i);
    a2 = sinGaus(a2');
    a3 = a2 * w3(:,:,i);
    a3 = mytanh(a3);
    ffdGenome(i,:) = a3(:)'; % Interleave vectors
end

foil = ffdRaeY(ffdGenome);
end

