clear all; close all; clc;
addpath('../DataScripts/SimulateData/InitFunctions')
addpath('../mdembedding')

%% Generate time series data

nvars = 20;
bc = 'fixed';

endtime = 10;
nobs = 500;
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.3;
cfn = @(n) constfn(n, damping);

prob = 0.5;
strength = 1;

mat = MakeNetworkSymmER(nvars, prob, true);
K = MakeNetworkTriDiag(nvars+2, false);
K(2:nvars+1, 2:nvars+1) = mat;
K = strength * K;

forcingFunc = zeros([nvars, length(tSpan)]);
%forcingFunc(ceil(nvars/2), 1) = 10;

numTrials = 1;
data = GenerateHarmonicData(nvars, tSpan, numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
noisyData = noisefn(data);


%% Use TISEAN to find optimal embedding dimension

tau1 = mdDelay(data.', 'maxLag', 50, 'plottype', 'all')
tau2 = mdDelay(data.', 'maxLag', 50, 'plottype', 'mean')
[fnnPercent, embeddingDimension] = mdFnn(data.', round(tau1))
