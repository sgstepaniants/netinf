clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'PertVarySizeForcingStrengths';

% Network sizes
networkSizes = 20;
numSizes = length(networkSizes);

% Connection strengths
strengths = 100; % IMPORTANT PARAMETER
numStrengths = length(strengths);

% Forcing magnitudes
forces = 100; % IMPORTANT PARAMETER
numForces = length(forces);

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Specify noise and prepocessing for data.
measParam = 0.01; % CONSTANT PARAMETER
noisefn = @(data) WhiteGaussianNoise(data, measParam);

% Delta t
deltat = 0.01; % CONSTANT PARAMETER

% Preprocessing function for data.
preprocfn = @(data) data;

% Amount of time to wait between perturbations
waitTime = 10;

% Probabilities of network connections.
prob = 0.5;

% Number of matrices to average results over.
numMats = 2;

% Number of experimental trials
numTrials = 1;

% Width of moving variance window
movvarWidth = 10;

% Threshold for mean variance algorithm
meanThresh = 0; %0.01; %0.001; % IMPORTANT PARAMETER

method = 'meanvar';

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);
if exist(expPath, 'dir') == 7
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
    if upper(m) == 'N'
        return
    end
    rmdir(expPath, 's')
end
mkdir(expPath)

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/PertResults', expPath);
if exist(resultPath, 'dir') == 7
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
    rmdir(resultPath, 's')
end
mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

% Run PCI to infer network connections.
predMats = cell(1, numSizes * numForces * numStrengths * numMats);
tprLog = nan(1, numSizes * numForces * numStrengths * numMats);
fprLog = nan(1, numSizes * numForces * numStrengths * numMats);
accLog = nan(1, numSizes * numForces * numStrengths * numMats);
numRerun = zeros(1, numSizes * numForces * numStrengths * numMats);

parsave = @(fname, noisyData, pertLength, pertTimes, mat)...
            save(fname, 'noisyData', 'pertLength', 'pertTimes', 'mat');

% Number of parallel processes
M = 25;
c = progress(numSizes * numForces * numStrengths * numMats);
parfor (idx = 1 : numSizes * numForces * numStrengths * numMats, M)
    [j, k, l, m] = ind2sub([numSizes, numForces, numStrengths, numMats], idx);
    fprintf('size: %d, force: %d, strength: %d\n', j, k, l)
    
    currExpPath = sprintf('%s/size%d/force%d/strength%d/mat%d', expPath, j, k, l, m);
    if exist(currExpPath, 'dir') ~= 7
        mkdir(currExpPath)
    end
    
    % Count the number of iterations done by the parfor loop
    c.count();

    nvars = networkSizes(j);
    force = forces(k);
    strength = strengths(l);

    % Create adjacency matrices.
    mat = MakeNetworkER(nvars, prob, true);

    % Perturb all nodes sequentially.
    pertIdx = 1 : nvars;
    numPerts = length(pertIdx);

    endtime = waitTime * (numPerts + 1);
    nobs = round(endtime / deltat);
    tSpan = linspace(0, endtime, nobs);

    % Build up forcing function.
    times = round(linspace(0, nobs, numPerts+2));
    pertTimes = times(2:end-1);
    pertLength = round(waitTime / (10 * deltat));

    forcingFunc = zeros([nvars, nobs]);
    for p=1:numPerts
        forcingFunc(pertIdx(p), pertTimes(p):pertTimes(p)+pertLength) = force;
    end

    % Generate data with forced perturbations.
    data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
    noisyData = noisefn(data);

    obsIdx = true([1, nvars]);
    leftPad = 0;
    rightPad = pertLength;
    truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
    [est, ~, ~, ~, ~, tableResults] = ...
        PerturbationBaseExperiment(noisyData, mat, numTrials, preprocfn, ...
            obsIdx, pertIdx, pertTimes, leftPad, rightPad, method, meanThresh, movvarWidth);

    parsave(sprintf('%s/dataLog.mat', currExpPath), noisyData, pertLength, pertTimes, mat);
    
    predMats{idx} = est;
    tprLog(idx) = tableResults.tpr;
    fprLog(idx) = tableResults.fpr;
    accLog(idx) = tableResults.acc;
end


predMats = reshape(predMats, numSizes, numForces, numStrengths, numMats);
tprLog = reshape(tprLog, [numSizes, numForces, numStrengths, numMats]);
fprLog = reshape(fprLog, [numSizes, numForces, numStrengths, numMats]);
accLog = reshape(accLog, [numSizes, numForces, numStrengths, numMats]);
numRerun = sum(reshape(numRerun, [numSizes, numForces, numStrengths, numMats]), 4);

% Save experiment results
save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/tprLog.mat', resultPath), 'tprLog');
save(sprintf('%s/fprLog.mat', resultPath), 'fprLog');
save(sprintf('%s/accLog.mat', resultPath), 'accLog');
save(sprintf('%s/numRerun.mat', resultPath), 'numRerun');


%% Plot Results

forceInd = 5;

% Show number of simulations that were skipped.
figure(1)
imagesc(reshape(numRerun(:, forceInd, :), [numSizes, numStrengths]))
set(gca,'YDir','normal')
colormap jet
colorbar
title('Number of Simulations Rerun by Our Analysis')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTick', strengths)
set(gca, 'YTick', networkSizes)
set(gca,'TickLength', [0 0])


% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 4);
figure(2)
clims = [0, 1];
imagesc(reshape(aveAccuracies(:, forceInd, :), [numSizes, numStrengths]), clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average Accuracy over Simulations')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTick', strengths)
set(gca, 'YTick', networkSizes)
%set(gca, 'TickLength', [0 0])


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 4);
figure(3)
clims = [0, 1];
imagesc(reshape(aveTPR(:, forceInd, :), [numSizes, numStrengths]), clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average TPR over Simulations')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTick', strengths)
set(gca, 'YTick', networkSizes)
%set(gca, 'TickLength', [0 0])


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 4);
figure(4)
clims = [0, 1];
imagesc(reshape(aveFPR(:, forceInd, :), [numSizes, numStrengths]), clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average FPR over Simulations')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTick', strengths)
set(gca, 'YTick', networkSizes)
%set(gca, 'TickLength', [0 0])
