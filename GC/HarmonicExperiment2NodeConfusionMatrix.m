clear all; close all; clc;
run('../mvgc_v1.0/startup.m')
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'ConfusionMatrix_Size2';

% Make directory to hold data files if one does not already exist
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../HarmonicExperiments/%s', expName);
if exist(expPath, 'dir') ~= 7
    mkdir(expPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', dataPath),'s');
    if upper(m) == 'N'
       return
    end
end

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/GCResults', expPath);
if exist(resultPath, 'dir') ~= 7
    mkdir(resultPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
end


%% Simulate Harmonic Oscillator Trajectories

% Network size
nvars = 2;

% Initialize masses, positions, and velocities of oscillators.
mfn = @(n) constfn(n, 1);
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) zeros([n, 1]);

% Specify the damping constant.
damping = 0;
cfn = @(n) constfn(n, damping);

% Define time sampling.
deltat = 0.1; % space between time points
endtime = 25;
nobs = round(endtime / deltat); % number of time points (observations)
tSpan = linspace(0, endtime, nobs);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

% Specify forcing function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Specify boundary conditions.
bc = 'fixed';

% Number of connectivity matrices to generate.
s = 100;
numMats = 4 * s;

% Number of simulation trials.
numTrials = 100;

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));

% Create 4 connectivity matrices for two oscillators.
trueMats = nan(nvars, nvars, numMats);
trueMats(:, :, 1:s)   = repmat([0, 0; 0, 0], [1, 1, s]);   % no connections
trueMats(:, :, s+1:2*s) = repmat([0, 0; 1, 0], [1, 1, s]);   % node 1 causes node 2
trueMats(:, :, 2*s+1:3*s) = repmat([0, 1; 0, 0], [1, 1, s]);   % node 2 causes node 1
trueMats(:, :, 3*s+1:4*s) = repmat([0, 1; 1, 0], [1, 1, s]);   % both nodes cause each other

% Turn these connectivity matrices into matrices of spring constants (connection strengths).
strengths = repmat(linspace(10/s, 10, s), [1, 4]);
Ks = nan(nvars+2, nvars+2, numMats);
for j = 1 : numMats
    K = MakeNetworkTriDiag(nvars+2, false);
    K(2:nvars+1, 2:nvars+1) = trueMats(:, :, j);
    K = strengths(j) * K;
    Ks(:, :, j) = K;
end

% Simulate oscillator trajectories.
dataLog = nan(nvars, nobs, numTrials, numMats);
for j = 1 : numMats
    fprintf('Computing simulations for matrix %d\n', j)
    K = Ks(:, :, j);
    data = GenerateHarmonicData(nvars, tSpan, ...
        numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
    noisyData = noisefn(data);
    dataLog(:, :, :, j) = noisyData;
end

% Save experiment simulated data and connectivity matrices.
save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/trueMats.mat', expPath), 'trueMats');
save(sprintf('%s/Ks.mat', expPath), 'Ks');


%% Run Granger Causality Experiments

% Run Granger Causality to infer network connections.
preprocfn = @(data) standardize(data);
save(sprintf('%s/expParams.mat', resultPath), 'preprocfn')

[predMats, tableResults] = GrangerBaseExperiment(dataLog, trueMats, preprocfn);

% Create a confusion matrix for network predictions.
confusionMat = ConfusionMatrix(trueMats, predMats)

save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/tableResults.mat', resultPath), 'tableResults');
save(sprintf('%s/confusionMat.mat', resultPath), 'confusionMat');
