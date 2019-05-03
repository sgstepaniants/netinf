% Simulates two connected spring mass oscillators with fixed boundary
% conditions. They have random initial positions and are released to
% oscillate.

clear all; close all; clc;
addpath('SimulateData/')
addpath('SimulateData/InitFunctions/')

expNum = 'Paper1';

% Run MVGC toolbox 10 times and vote for final connectivity matrix.
nvars = 2;

% Initialize masses, positions, and velocities of oscillators.
pfn = @(n) randfn(n, -0.5, 0.5);
wfn = @(n) randfn(n, -1, 1);

% Define time sampling.
deltat = 0.1; % space between time points
endtime = 25;
nobs = round(endtime / deltat); % number of time points (observations)
tSpan = linspace(0, endtime, nobs);

% Specify noise and prepocessing for data.
measParam = 0.01;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

% Specify forcing function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Number of connectivity matrices to generate.
numMats = 400;

% Number of simulation trials.
numTrials = 100;

% Make directory to hold data files if one does not already exist
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);
if exist(expPath, 'dir') ~= 7
    mkdir(expPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
    if upper(m) == 'N'
       return
    end
end

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));

% Create 4 connectivity matrices for two oscillators.
trueMats = nan(nvars, nvars, numMats);
trueMats(:, :, 1:100)   = repmat([0, 0; 0, 0], [1, 1, 100]);   % no connections
trueMats(:, :, 101:200) = repmat([0, 0; 1, 0], [1, 1, 100]);   % node 1 causes node 2
trueMats(:, :, 201:300) = repmat([0, 1; 0, 0], [1, 1, 100]);   % node 2 causes node 1
trueMats(:, :, 301:400) = repmat([0, 1; 1, 0], [1, 1, 100]);   % both nodes cause each other

% Make a vector of connection strengths.
Ks = repmat(0.1 : 0.1 : 10, [1, 4]);

% Simulate oscillator trajectories.
dataLog = nan(nvars, nobs, numTrials, numMats);
for j = 1 : numMats
    fprintf('Computing simulations for matrix %d\n', j)
    A = trueMats(:, :, j);
    K = Ks(j);
    data = GenerateKuramotoData(A, tSpan, ...
        numTrials, K, pfn, wfn, forcingFunc);
    noisyData = noisefn(data);
    dataLog(:, :, :, j) = noisyData;
end

% Save experiment simulated data and connectivity matrices.
save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/trueMats.mat', expPath), 'trueMats');
save(sprintf('%s/Ks.mat', expPath), 'Ks');
