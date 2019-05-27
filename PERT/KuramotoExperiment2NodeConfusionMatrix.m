clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')
addpath('../BAgraph_dir/')
addpath('../kmeans_opt/')

%% Initialize Parameters

% Number of network nodes
nvars = 2;

s = 100;
numMats = 4 * s;

% Create 4 connectivity matrices for two oscillators.
mats = zeros(nvars, nvars, numMats);
mats(:, :, 1:s)   = repmat([0, 0; 0, 0], [1, 1, s]);   % no connections
mats(:, :, s+1:2*s) = repmat([0, 0; 1, 0], [1, 1, s]);   % node 1 causes node 2
mats(:, :, 2*s+1:3*s) = repmat([0, 1; 0, 0], [1, 1, s]);   % node 2 causes node 1
mats(:, :, 3*s+1:4*s) = repmat([0, 1; 1, 0], [1, 1, s]);   % both nodes cause each other
numMats = size(mats, 3);

% Turn these connectivity matrices into matrices of connection strengths.
strengths = repmat(linspace(20/s, 20, s), [1, 4]) + 20;

% Gaussian noise function
noiseVar = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, noiseVar);

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Delta t
deltat = 0.1;

% Perturbation force for oscillators
pertForce = 30;

% Amount of time to wait between perturbations
waitTime = 3;

% Width of moving variance window
movvarWidth = round(waitTime / (20 * deltat));

% Threshold for mean variance algorithm
meanThresh = 0.01;

% Padding for window
pad = 0;

% Number of experimental trials
numTrials = 1;

method = 'meanvar';

% Make directory to hold results files if one does not already exist
expName = 'PertConfusionMatrix_Size2';
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

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/PertResults', expPath);
if exist(resultPath, 'dir') ~= 7
    mkdir(resultPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
end


%% Simulate Data with Varying Connectivity Matrices
dataLog = cell(1, numMats);
dataPertIdx = zeros(numMats, nvars);
dataPertTimes = zeros(numMats, nvars);
dataPertLength = zeros(1, numMats);

fprintf('Simulate Data:\n')
for j = 1 : numMats
    j
    % Perturb all nodes sequentially.
    pertIdx = 1:nvars;
    numPerts = length(pertIdx);

    % Build up network connectivity
    mat = mats(:, :, j);
    strength = strengths(j);
    
    % Create the timespan for the simulation.
    endtime = waitTime * (numPerts + 1);
    nobs = round(endtime / deltat);
    tSpan = linspace(0, endtime, nobs);

    % Build up forcing function.
    times = round(linspace(0, nobs, numPerts+2));
    pertTimes = times(2:end-1);
    pertLength = round(nobs/(10*(numPerts+1)));

    forcingFunc = zeros([nvars, nobs]);
    for k=1:numPerts
        forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
    end

    % Generate data with forced perturbations.
    data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
    noisyData = noisefn(data);

    dataLog{j} = noisyData;
    dataPertIdx(j, :) = pertIdx;
    dataPertTimes(j, :) = pertTimes;
    dataPertLength(j, :) = pertLength;

    % Save experiment simulated data, connectivity matrices, and parameters.
    save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
    save(sprintf('%s/dataPertIdx.mat', expPath), 'dataPertIdx');
    save(sprintf('%s/dataPertTimes.mat', expPath), 'dataPertTimes');
    save(sprintf('%s/dataPertLength.mat', expPath), 'dataPertLength');
end


%% Evaluate Algorithm on Data and Create Confusion Matrix

% Nodes to observe.
dataObsIdx = true(numMats, 2);
preprocfn = @(data) standardize(data);

fprintf('Run Algorithm:\n')
[predMats, predMatsHist, AprobHist, truePertOrders, predPertOrders] = ...
    PerturbationBaseExperiment(dataLog, mats, numTrials, preprocfn, ...
        dataObsIdx, dataPertIdx, dataPertTimes, dataPertLength, method, ...
        meanThresh, pad, movvarWidth);

predMats1 = predMatsHist(:, :, 1, :);
predMats1(isnan(predMats1)) = 0;
predMats2 = predMatsHist(:, :, 2, :);
predMats2(isnan(predMats2)) = 0;

confusionMatPert1 = ConfusionMatrix(mats, predMats1)
confusionMatPert2 = ConfusionMatrix(mats, predMats2)

save(sprintf('%s/confusionMatPert1.mat', resultPath), 'confusionMatPert1');
save(sprintf('%s/confusionMatPert2.mat', resultPath), 'confusionMatPert2');
save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/predMatsHist.mat', resultPath), 'predMatsHist');
save(sprintf('%s/AprobHist.mat', resultPath), 'AprobHist');
