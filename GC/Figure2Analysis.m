clear all; close all; clc;
addpath('../DataScripts/SimulateData/InitFunctions')
addpath('../mdembedding')

%% Generate time series data

% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

nvars = 5;

endtime = 5;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

prob = 0.5;
%mat = zeros(nvars);
%mat(1:(nvars/2), 1:(nvars/2)) = ones(nvars/2);
%mat((nvars/2)+1:nvars, (nvars/2)+1:nvars) = ones(nvars/2);
%mat = mat - eye(nvars);
mat = MakeNetworkER(nvars, prob, true);

strength = 5;

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Specify noise for data
measParam = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, measParam);

forcingFunc = zeros(nvars, nobs);

numTrials = 100;
data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
noisyData = noisefn(data);

figure(1)
plot(noisyData(:, :, 1).')


%% Run Granger Causality on Data
rhoThresh = 1;
[predMat, diagnostics] = DemoMVGC(cos(noisyData), rhoThresh);
figure(2)
imagesc(mat)
figure(3)
imagesc(predMat)
acc = 1 - nnz(predMat - mat) / (nvars^2 - nvars)


%% Create Random Matrices with Certain Predictive Accuracy (For Instructive Purposes Only)
inds = 1 : nvars^2;
inds((nvars + 1) * (1 : nvars) - nvars) = [];
maxEdges = nvars^2 - nvars;

refAcc = 0.9;
% Edges that we will change randomly and uniformly in the original network
% to get the prescribed accuracy.
wrongInds = randsample(inds, round((1 - acc) * maxEdges));

predMatRand90 = mat;
predMatRand90(wrongInds) = ~predMatRand90(wrongInds)
