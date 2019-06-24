clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions')
addpath('../mdembedding')

%% Generate time series data

% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

nvars = 8;

expNum = 'NetworkInferenceAccuracy';

prob = 0.5;

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Specify noise and preprocessing function for data
preprocfn = @(data) cos(data);
forcingFunc = zeros(nvars, nobs);

numTrials = 100;

rhoThresh = 1;

% Check that directory with experiment data exists
%expName = sprintf('EXP%s', expNum);
%expPath = sprintf('../KuramotoExperiments/%s', expName);
%if exist(expPath, 'dir') == 7
%    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
%    if upper(m) == 'N'
%        return
%    end
%    rmdir(expPath, 's')
%end
%mkdir(expPath)

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));

mat = zeros(nvars);
mat(1:(nvars/2), 1:(nvars/2)) = ones(nvars/2);
mat((nvars/2)+1:nvars, (nvars/2)+1:nvars) = ones(nvars/2);
mat = mat - eye(nvars);
matEigs = eig(mat);

% Save matrix
save(sprintf('%s/trueMat.mat', expPath), 'mat');

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/GCResults', expPath);
if exist(resultPath, 'dir') == 7
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
    rmdir(resultPath, 's')
end
mkdir(resultPath)


%% Run Granger Causality on Data

load(sprintf('%s/dataLog.mat', expPath))
load(sprintf('%s/predMats.mat', resultPath))

% Granger network with 90% accuracy
%endtime = 5;
%deltat = 0.1;
%nobs = round(endtime / deltat);
%tSpan = linspace(0, endtime, nobs);
%forcingFunc = zeros(nvars, nobs);
%data1 = GenerateKuramotoData(mat, tSpan, numTrials, 8, pfn, wfn, forcingFunc);
%noisyData1 = WhiteGaussianNoise(data1, 0.1);

%[predMat1, diagnostics1] = DemoMVGC(preprocfn(noisyData1), rhoThresh);
predMat1Eigs = eig(predMat1);
acc1 = 1 - nnz(predMat1 - mat) / (nvars^2 - nvars)


% Granger network with 80% accuracy
%endtime = 5;
%deltat = 0.1;
%nobs = round(endtime / deltat);
%tSpan = linspace(0, endtime, nobs);
%forcingFunc = zeros(nvars, nobs);
%data2 = GenerateKuramotoData(mat, tSpan, numTrials, 2, pfn, wfn, forcingFunc);
%noisyData2 = WhiteGaussianNoise(data2, 0.1);

%[predMat2, diagnostics2] = DemoMVGC(preprocfn(noisyData2), rhoThresh);
predMat2Eigs = eig(predMat2);
acc2 = 1 - nnz(predMat2 - mat) / (nvars^2 - nvars)


% Granger network with 75% accuracy
%endtime = 25;
%deltat = 0.1;
%nobs = round(endtime / deltat);
%tSpan = linspace(0, endtime, nobs);
%forcingFunc = zeros(nvars, nobs);
%data3 = GenerateKuramotoData(mat, tSpan, numTrials, 8, pfn, wfn, forcingFunc);
%noisyData3 = WhiteGaussianNoise(data3, 0.1);

%[predMat3, diagnostics3] = DemoMVGC(preprocfn(noisyData3), rhoThresh);
predMat3Eigs = eig(predMat3);
acc3 = 1 - nnz(predMat3 - mat) / (nvars^2 - nvars)

%save(sprintf('%s/dataLog.mat', expPath), 'noisyData1', 'noisyData2', 'noisyData3')
save(sprintf('%s/predMats.mat', resultPath), 'predMat1', 'predMat2', 'predMat3')

figure(1)
plot(real(matEigs), imag(matEigs), 'b*')
figure(2)
plot(real(predMat1Eigs), imag(predMat1Eigs), 'g*')
figure(3)
plot(real(predMat2Eigs), imag(predMat2Eigs), 'y*')
figure(4)
plot(real(predMat3Eigs), imag(predMat3Eigs), 'r*')

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
