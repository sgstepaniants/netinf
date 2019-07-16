% Generate time series data for n masses connected to their nearest
% neighbors by springs with certain initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')
addpath('../BAgraph_dir/')

warning('off', 'stats:kmeans:EmptyCluster')

nvars = 10;
bc = 'fixed';

endtime = 200;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn  = @(data) WhiteGaussianNoise(data, 0.1);

pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.2;
cfn = @(n) constfn(n, damping);

% Create connectivity matrix.
prob = 0.1;
spring = 0.1;
%mat = MakeNetworkSymmER(nvars, prob, true);
%mat(2:end, 1) = 1;
mat = BAgraph_dir(nvars, floor(nvars/2), floor(nvars/2));
K = MakeNetworkTriDiag(nvars + 2, false);
K(2:nvars+1, 2:nvars+1) = mat;
K = spring * K;

% Create forcing function.
forcingFunc = zeros([nvars, length(tSpan)]);
pertIdx = 8;
numPerts = length(pertIdx);
times = round(linspace(0, length(tSpan), numPerts+2));
pertTimes = times(2:end-1);
pertLength = 100;
pertForce = 10;
for k=1:numPerts
    forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
end

% Simulate harmonic oscillator movement.
ntrials = 1;
data = GenerateHarmonicData(nvars, tSpan, 1, K, pfn, vfn, mfn, cfn, bc, forcingFunc);

% Plot oscillator trajectories.
figure(1)
plot(data.');
legend(strcat('n', num2str((1:nvars).')))
