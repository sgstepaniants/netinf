clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions')
addpath('../mdembedding')

%% Generate time series data

% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

nvars = 5;

% Initial conditions
pfn = @(n) [0; -1; -1; 1; 2];
wfn = @(n) linspace(0, 2*pi, n).';

K = 7;

mat1 = zeros(nvars);
mat1(2, 1) = 1;
mat1(3, 1) = 1;
mat1(4, 2) = 1;
mat1(2, 3) = 1;
mat1(5, 3) = 1;

mat2 = zeros(nvars);
mat2(2, 1) = 1;
mat2(3, 1) = 1;
mat1(2, 3) = 1;
mat2(4, 3) = 1;
mat2(5, 3) = 1;

endtime = 10;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);
forcingFunc = zeros(nvars, nobs);

data1 = GenerateKuramotoData(mat1, tSpan, 1, K, pfn, wfn, forcingFunc);
figure(1)
plot(data1.')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'TickLength', [0 0])

data2 = GenerateKuramotoData(mat2, tSpan, 1, K, pfn, wfn, forcingFunc);
figure(2)
plot(data2.')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'TickLength', [0 0])
