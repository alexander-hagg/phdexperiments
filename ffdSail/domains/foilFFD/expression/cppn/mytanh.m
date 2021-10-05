function [ y ] = mytanh( x )
%MYTANH Custon sigmoid activation function
%   y = 1./(1+exp(-2*(x*2.5)))

y = 1./(1+exp(-2*(x*2.5)));


end

