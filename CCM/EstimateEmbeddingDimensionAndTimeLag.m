clear all; close all; clc;
addpath('../DataScripts/SimulateData')
addpath('../DataScripts/SimulateData/InitFunctions')
addpath('../mdembedding')

%% Generate time series data

nvars = 30;
bc = 'fixed';

maxDelay = 100;
maxEmb = 20;

endtime = 25;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0;
cfn = @(n) constfn(n, damping);

prob = 0.1;
strength = 1;

mat = MakeNetworkSymmER(nvars, prob, true);
K = MakeNetworkTriDiag(nvars+2, false);
K(2:nvars+1, 2:nvars+1) = mat;
K = strength * K;

forcingFunc = zeros([nvars, length(tSpan)]);

numTrials = 1;
data = GenerateHarmonicData(nvars, tSpan, numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);

noisefn  = @(data) WhiteGaussianNoise(data, 0.1);
noisyData = noisefn(data);


%% Use mdembedd to find optimal time lag and embedding dimension

tau = round(mdDelay(data.', 'maxLag', maxDelay, 'plottype', 'mean'))
currMaxEmb = nobs - 1;
if tau > 1; currMaxEmb = floor(nobs / tau); end
[fnnPercent, Es] = mdFnn(data(1, :).', tau, 'maxEmb', min(currMaxEmb, maxEmb), 'doPlot', 1);
E = findElbow(Es, fnnPercent)
