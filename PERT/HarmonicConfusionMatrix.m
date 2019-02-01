clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')
addpath('./BAgraph_dir/')
addpath('./kmeans_opt/')

%% Initialize Parameters

% Number of network nodes
nvars = 2;
% Boundary conditions of oscillator system
bc = 'fixed';

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

% Gaussian noise function
noiseVar = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, noiseVar);

% Spring constants in oscillator network
springs = 0.1:0.1:1;
% Damping in oscillator network
damping = 0.3;

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);
cfn = @(n) constfn(n, damping);

% Number of nodes we can observe in the network
numObs = nvars;

% Delta t
deltat = 0.1;

% Perturbation force for oscillators
pertForce = 30;

% Threshold for correlation algorithm
corrThresh = 0.5;

% Padding for window
pad = 100;

% Number of experimental trials
numTrials = 100;

method = 'corr';

% Make directory to hold results files if one does not already exist
expName = sprintf('EXPConfusionMatrix(nvars%d_numObs%d_damping%.2f_pertf%.2f)', nvars, numObs, damping, pertForce);
expPath = sprintf('HarmonicExperiments/%s', expName);
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


%% Siumlate Data with Varying Connectivity Matrices
dataLog = cell(1, numMats);
dataParamsLog = cell(numMats, 3);

fprintf('Simulate Data:\n')
for j = 1:numMats
    j
    % Perturb all nodes sequentially.
    pertIdx = 1:nvars;
    numPerts = length(pertIdx);

    % Build up network connectivity
    mat = mats(:, :, j);
    K = Ks(:, :, j);

    % Check if this adjacency matrix has disconnected oscillators.
    G = digraph(mat.');
    distLeft = distances(G, 1);
    distRight = distances(G, nvars);
    disconnectedNodes = find(~isfinite(distLeft) & ~isfinite(distRight));

    % Check if this adjacency matrix is resonant.
    A = K(2:nvars+1, 2:nvars+1);
    A = A - diag(sum(K(2:nvars+1, :), 2));
    lambdas = eig(A);
    amplitudes = real(-damping + sqrt(damping^2 + 4 * lambdas)) / 2;

    % If this adjacency matrix is bad, make a new simulation.
    if ~isempty(disconnectedNodes) || any(amplitudes > -0.00001)
        continue
    end

    % Create the timespan for the simulation.
    eps = 0.01;
    waitTime = ceil(log(eps / min(sqrt(sum((pertForce * inv(A)).^2)))) / max(amplitudes));
    if waitTime > 500
        continue
    end

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
    data = GenerateNNCoupledData(nvars, tSpan, 1, K, pfn, vfn, ...
        mfn, cfn, bc, forcingFunc);
    noisyData = noisefn(data);

    dataLog{j} = noisyData;
    dataParamsLog{j, 1} = pertIdx;
    dataParamsLog{j, 2} = pertTimes;
    dataParamsLog{j, 3} = pertLength;

    % Save experiment simulated data, connectivity matrices, and parameters.
    save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
    save(sprintf('%s/dataParamsLog.mat', expPath), 'dataParamsLog');
end


%% Evaluate Algorithm on Data and Create Confusion Matrix
truePertOrdersLog = zeros(nvars, nvars, numMats);
predPertOrdersLog = zeros(nvars, nvars, numMats);
predMatsHist = zeros(nvars, nvars, numMats, nvars);

% Nodes to observe.
obsIdx = logical([1, 1]);

fprintf('Run Algorithm:\n')
for j = 1:numMats
    data = dataLog{j};
    mat = mats(:, :, j);
    pertIdx = dataParamsLog{j, 1};
    pertTimes = dataParamsLog{j, 2};
    pertLength = dataParamsLog{j, 3};

    % Select only the subset of nodes that we can observe.
    observedData = data(obsIdx, :);
    [AprobHist, predPertOrders] = CreateProbabilityMatrix(observedData, ...
                                                    pertIdx, obsIdx, pertTimes, ...
                                                    pertLength, method, corrThresh, pad);
    predPertOrdersLog(:, :, j) = predPertOrders;

    % Get the true perturbation orders.
    truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
    truePertOrdersLog(:, :, j) = truePertOrders;

    % Get the network reconstruction our algorithm produces.
    predMatHist = AprobHist > 0.5;
    predMatsHist(:, :, j, :) = predMatHist;

    % Save experiment results.
    save(sprintf('%s/truePertOrdersLog.mat', expPath), 'truePertOrdersLog');
    save(sprintf('%s/predPertOrdersLog.mat', expPath), 'predPertOrdersLog');
    save(sprintf('%s/predMatsHist.mat', expPath), 'predMatsHist');
end

confusionMatPert1 = ConfusionMatrix(mats, predMatsHist(:, :, :, 1))
confusionMatPert2 = ConfusionMatrix(mats, predMatsHist(:, :, :, 2))
