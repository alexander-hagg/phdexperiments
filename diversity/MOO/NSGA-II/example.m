% ----------------------------------------------------------------------- %
% Example of use of the funcion NSGAII.m, which performs a Non Sorting Gene-
% tic Algorithm II (NSGA-II), based on Deb2002.                           %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez Cagigal                                      %
%   Date:    22/12/2017                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%    [1] Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. A. M. T. (2002)%
%        A fast and elitist multiobjective genetic algorithm: NSGA-II.    %
%        IEEE transactions on evolutionary computation, 6(2), 182-197.    %
% ----------------------------------------------------------------------- %
clear all; clc; close all;


MultiObj.fun = @(x) [-10.*(exp(-0.2.*sqrt(x(:,1).^2+x(:,2).^2)) + exp(-0.2.*sqrt(x(:,2).^2+x(:,3).^2))), ...
    sum(abs(x).^0.8 + 5.*sin(x.^3),2)];
MultiObj.nVar = 3;
MultiObj.var_min = -5.*ones(1,MultiObj.nVar);
MultiObj.var_max = 5.*ones(1,MultiObj.nVar);
%load('Kursawe.mat');
%MultiObj.truePF = PF;

% Parameters
params.Np = 10*10;        % Population size
params.pc = 0.9;        % Probability of crossover
params.pm = 0.5;        % Probability of mutation
params.maxgen = 500;    % Maximum number of generations

% MOPSO
profile on
[R, Rfit] = NSGAII(params,MultiObj);
profile off
profile viewer

