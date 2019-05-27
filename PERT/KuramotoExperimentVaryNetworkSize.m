clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'PertVaryNetworkSize';

networkSizes = 5;
numSizes = length(networkSizes);

% Initialize masses, positions, and velocities of oscillators.
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Define time sampling.
deltat = 0.1; % space between time points

% Specify noise and prepocessing for data.
noiseVar = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, noiseVar);

% Perturbation force for oscillators
pertForce = 50;

% Amount of time to wait between perturbations
waitTime = 3;

% Probabilities of network connections.
prob = 0.5;

% Connection strength.
strength = 50;

% Number of matrices to average results over.
numMats = 10;

% Number of simulation trials per repetition.
numTrials = 1;

% Width of moving variance window
movvarWidth = round(waitTime / (20 * deltat));

% Padding for window
pad = 0;

% Threshold for mean variance algorithm
meanThresh = 0.01;

method = 'meanvar';

% Check that directory with experiment data exists
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

%% Generate Data and Run Granger Causality Experiments

% Create random connectivity matrices and simulate oscillator trajectories.
dataLog = cell(1, numSizes * numMats);
trueMats = cell(1, numSizes * numMats);

% Run PCI to infer network connections.
preprocfn = @(data) cos(data);
save(sprintf('%s/expParams.mat', resultPath), 'preprocfn')

predMats = cell(1, numSizes * numMats);
tprLog = nan(1, numSizes * numMats);
fprLog = nan(1, numSizes * numMats);
accuracyLog = nan(1, numSizes * numMats);

% Number of parallel processes
M = 12;
c = progress(numSizes * numMats);
parfor (idx = 1 : numSizes * numMats, M)
    [j, l] = ind2sub([numSizes, numMats], idx);
    fprintf('size: %d, mat: %d\n', j, l)
    
    % Count the number of iterations done by the parfor loop
    c.count();
    
    nvars = networkSizes(j);
    
    % Create adjacency matrices.
    mat = MakeNetworkER(nvars, prob, true);
        
    % Perturb all nodes sequentially.
    pertIdx = 1:nvars;
    numPerts = length(pertIdx);
    
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

    obsIdx = true([1, nvars]); % default parameter
    [est, predMatsHist, AprobHist, truePertOrders, predPertOrders, tableResults] = ...
        PerturbationBaseExperiment(data, mat, numTrials, preprocfn, ...
            obsIdx, pertIdx, pertTimes, pertLength, method, ...
            meanThresh, pad, movvarWidth);
    
    dataLog{idx} = noisyData;
    trueMats{idx} = mat;
    predMats{idx} = est;

    tprLog(idx) = tableResults.tpr;
    fprLog(idx) = tableResults.fpr;
    accuracyLog(idx) = tableResults.acc;
end

% Reshape data structures
dataLog = reshape(dataLog, numSizes, numMats);
trueMats = reshape(trueMats, numSizes, numMats);
predMats = reshape(predMats, numSizes, numMats);
tprLog = reshape(tprLog, [numSizes, numMats]);
fprLog = reshape(fprLog, [numSizes, numMats]);
accuracyLog = reshape(accuracyLog, [numSizes, numMats]);

% Save experiment simulated data and connectivity matrices.
save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/trueMats.mat', expPath), 'trueMats');

% Save experiment results
save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/tprLog.mat', resultPath), 'tprLog');
save(sprintf('%s/fprLog.mat', resultPath), 'fprLog');
save(sprintf('%s/accuracyLog.mat', resultPath), 'accuracyLog');


%% Plot Results
% 
% % Show average accuracies for each number of perturbations and
% % observations.
% aveAccuracies = nanmean(accuracyLog, 2);
% figure(2)
% plot(networkSizes, aveAccuracies)
% title('Average Accuracy over Simulations')
% xlabel('Network Size')
% ylabel('Accuracy')
% xlim([networkSizes(1) networkSizes(end)]) 
% ylim([0 1])
% 
% % Show average TPR for each number of perturbations and
% % observations.
% aveTPR = nanmean(tprLog, 2);
% figure(3)
% plot(networkSizes, aveTPR)
% title('Average TPR over Simulations')
% xlabel('Network Size')
% ylabel('TPR')
% xlim([networkSizes(1) networkSizes(end)]) 
% ylim([0 1])
% 
% % Show average FPR for each number of perturbations and
% % observations.
% aveFPR = nanmean(fprLog, 2);
% figure(4)
% plot(networkSizes, aveFPR)
% title('Average FPR over Simulations')
% xlabel('Network Size')
% ylabel('FPR')
% xlim([networkSizes(1) networkSizes(end)]) 
% ylim([0 1])
