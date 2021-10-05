function [ y ] = sinGaus( x )
%SINGAUS Custom activation function for sin of 1st act gaus of 2nd
% y(:,1) = 0.5+0.5*(sin(x(1,:)*2*pi));          % Custom Sine
% y(:,2) = (exp(-((x(2,:)*2.5)-0).^2/(2*1^2))); % Custom Gaus

y(:,1) = 0.5+0.5*(sin(x(1,:)*2*pi));
y(:,2) = (exp(-((x(2,:)*2.5)-0).^2/(2*1^2)));


end
