clear all; close all; clc;
addpath('../DataScripts/SimulateData/InitFunctions')
addpath('../mdembedding')

%% Generate time series data

% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

clear all; close all; clc;
addpath('./InitFunctions/')

nvars = 10;

endtime = 10;
nobs = 250;
tSpan = linspace(0, endtime, nobs);

prob = 0.5;
mat = zeros(nvars);
mat(1, 2) = 1;
mat(2, 3) = 1;
mat(4, 3) = 1;
mat(5, 3) = 1;
mat(4, 5) = 1;
mat(6, 7) = 1;
mat(6, 10) = 1;
mat(7, 10) = 1;
mat(8, 10) = 1;
mat(9, 10) = 1;

strength = 10;

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Specify noise for data
measParam = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, measParam);

forcingFunc = zeros(nvars, nobs);

numTrials = 1;
data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
noisyData = noisefn(data);

plot(noisyData.')


%% Run Granger Causality on Data
rhoThresh = 1;
[predMat, diagnostics] = DemoMVGC(cos(noisyData), rhoThresh);
