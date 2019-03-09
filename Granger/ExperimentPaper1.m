% Simulates two connected spring mass oscillators with fixed boundary
% conditions. They have random initial positions and are released to
% oscillate.

clear all; close all; clc;
addpath('../SimulateData/')

expnum = 'Paper1';

% Run MVGC toolbox 10 times and vote for final connectivity matrix.
nvars = 2;

% Initialize masses, positions, and velocities of oscillators.
mfn = @(n) constmfn(n, 1);
pfn = @(n) randpfn(n);
vfn = @(n) zeros([n, 1]);

% Specify the damping constant.
damping = 0.5;
cfn = @(n) constcfn(n, damping);

% Create 4 connectivity matrices for two oscillators.
mats = zeros(nvars, nvars, 400);
mats(:, :, 1:100)   = repmat([0, 0; 0, 0], [1, 1, 100]);   % no connections
mats(:, :, 101:200) = repmat([0, 0; 1, 0], [1, 1, 100]);   % node 1 causes node 2
mats(:, :, 201:300) = repmat([0, 1; 0, 0], [1, 1, 100]);   % node 2 causes node 1
mats(:, :, 301:400) = repmat([0, 1; 1, 0], [1, 1, 100]);   % both nodes cause each other
numMats = size(mats, 3);

% Turn these connectivity matrices into matrices of spring constants (connection strengths).
springConsts = repmat(0.1 : 0.1 : 10, [1, 4]);
Ks = zeros(nvars+2, nvars+2, numMats);
for j = 1 : numMats
    K = MakeNetworkTriDiag(nvars+2, false);
    K(2:nvars+1, 2:nvars+1) = mats(:, :, j);
    K = springConsts(j) * K;
    Ks(:, :, j) = K;
end

% Define time sampling.
deltat = 0.1; % space between time points
endtime = 25;
nobs = round(endtime / deltat); % number of time points (observations)
tSpan = linspace(0, endtime, nobs);

% Number of times to split time interval.
tsplits = nobs;
% How often to save worspace (after runs of MVGC).
freq = 1;

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);
preprocfn = @(data) NoiseThenDetrend(data, noisefn);

% Specify focring function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Specify boundary conditions.
bc = 'fixed';

% Number of simulation trials.
ntrials = 100;
% Number of times to repeat MVGC toolbox runs.
reps = 1;

% Simulate oscillator trajectories.
data = zeros(nvars, nobs, ntrials, reps, numMats);
for j = 1 : numMats
    fprintf('Computing simulations for matrix %d\n', j)
    K = Ks(:, :, j);
    for r = 1 : reps
        data(:, :, :, r, j) = GenerateNNCoupledData(nvars, tSpan, ...
            ntrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
    end
end

% Add noise and detrend oscillator data.
preprocessedData = preprocfn(data);

% Run Granger Causality to infer network connections.
BaseExperiment(expnum, mats, preprocessedData, ntrials, reps, tsplits, freq)

% Create a confusion matrix for network predictions.
load(sprintf('./exp%s/exp%s.mat', expnum, expnum))
predMats = squeeze(est);
confusionMat = ConfusionMatrix(mats, predMats)

save(sprintf('./exp%s/exp%s-params.mat', expnum, expnum))
