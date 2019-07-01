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

prob = 0.5;
strength = 10;

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
currMaxEmb = min(floor(nobs / tau), maxEmb);
[fnnPercent, Es] = mdFnn(data(1, :).', tau, 'maxEmb', currMaxEmb, 'doPlot', 1);

numEs = length(Es);
p1 = [Es(1), fnnPercent(1), 0];
p2 = [Es(end), fnnPercent(end), 0];
proj = cross(repmat(p2 - p1, [numEs, 1]), [Es.', fnnPercent.', zeros(numEs, 1)] - repmat(p1, [numEs, 1]), 2) / norm(p2 - p1);
dists = sqrt(sum(proj.^2, 2));
[~, E] = max(dists)
