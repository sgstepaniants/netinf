clear all; close all; clc;
addpath('../DataScripts/SimulateData')
addpath('../DataScripts/SimulateData/InitFunctions')
addpath('../mdembedding')

%% Generate time series data

nvars = 20;
bc = 'fixed';

endtime = 25;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.3);

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0;
cfn = @(n) constfn(n, damping);

prob = 0.5;
strength = 50;

mat = MakeNetworkSymmER(nvars, prob, true);
K = MakeNetworkTriDiag(nvars+2, false);
K(2:nvars+1, 2:nvars+1) = mat;
K = strength * K;

forcingFunc = zeros([nvars, length(tSpan)]);

numTrials = 1;
data = GenerateHarmonicData(nvars, tSpan, numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
noisyData = noisefn(data);


%% Use mdembedd to find optimal time lag and embedding dimension

tau = mdDelay(data.', 'maxLag', 50, 'plottype', 'mean')
fnnPercent = mdFnn(data(1, :).', round(tau), 'maxEmb', 10, 'doPlot', 1);
E = find(fnnPercent < 1, 1, 'first')
