clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'PertVarySizeForcingStrengths';

% Network sizes
networkSizes = 2 : 20;
numSizes = length(networkSizes);

% Connection strengths
strengths = 10 : 10 : 100; % IMPORTANT PARAMETER
numStrengths = length(strengths);

% Forcing magnitudes
forces = 10 : 10 : 100; % IMPORTANT PARAMETER
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
numMats = 100;

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
dataLog = cell(1, numSizes * numForces * numStrengths * numMats);
dataPertLength = cell(1, numSizes * numForces * numStrengths * numMats);
dataPertTimes = cell(1, numSizes * numForces * numStrengths * numMats);
trueMats = cell(1, numSizes * numForces * numStrengths * numMats);

% Run PCI to infer network connections.
predMats = cell(1, numSizes * numForces * numStrengths * numMats);
tprLog = nan(1, numSizes * numForces * numStrengths * numMats);
fprLog = nan(1, numSizes * numForces * numStrengths * numMats);
accLog = nan(1, numSizes * numForces * numStrengths * numMats);
numRerun = zeros(1, numSizes * numForces * numStrengths * numMats);

% Number of parallel processes
M = 25;
c = progress(numSizes * numForces * numStrengths * numMats);
for idx = 1 : numSizes * numForces * numStrengths * numMats %parfor (idx = 1 : numSizes * numForces * numStrengths * numMats, M)
    [j, k, l, m] = ind2sub([numSizes, numForces, numStrengths, numMats], idx);
    fprintf('size: %d, force: %d, strength: %d\n', j, k, l)
    
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
    for k=1:numPerts
        forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = force;
    end

    % Generate data with forced perturbations.
    data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
    noisyData = noisefn(data);

    obsIdx = true([1, nvars]);
    leftPad = 0;
    rightPad = pertLength;
    truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
    [est, ~, ~, ~, predPertOrders, tableResults] = ...
        PerturbationBaseExperiment(noisyData, mat, numTrials, preprocfn, ...
            obsIdx, pertIdx, pertTimes, leftPad, rightPad, method, meanThresh, movvarWidth);

    dataLog{idx} = noisyData;
    dataPertLength{idx} = pertLength;
    dataPertTimes{idx} = pertTimes;
    trueMats{idx} = mat;

    predMats{idx} = est;
    tprLog(idx) = tableResults.tpr;
    fprLog(idx) = tableResults.fpr;
    accLog(idx) = tableResults.acc;
end

% Reshape data structures
dataLog = reshape(dataLog, numSizes, numForces, numStrengths, numMats);
dataPertLength = reshape(dataPertLength, numSizes, numForces, numStrengths, numMats);
dataPertTimes = reshape(dataPertLength, numSizes, numForces, numStrengths, numMats);
trueMats = reshape(trueMats, numSizes, numForces, numStrengths, numMats);

predMats = reshape(predMats, numSizes, numForces, numStrengths, numMats);
tprLog = reshape(tprLog, [numSizes, numForces, numStrengths, numMats]);
fprLog = reshape(fprLog, [numSizes, numForces, numStrengths, numMats]);
accLog = reshape(accLog, [numSizes, numForces, numStrengths, numMats]);
numRerun = sum(reshape(numRerun, [numSizes, numForces, numStrengths, numMats]), 4);


% Save experiment data
for j = 1 : numSizes
    for k = 1 : numForces
        currExpPath = sprintf('%s/size%d/force%d', expPath, j, k);
        if exist(currExpPath, 'dir') ~= 7
            mkdir(currExpPath)
        end
        
        currDataLog = dataLog{j, k, :, :};
        currDataPertLength = dataPertLength{j, k, :, :};
        currDataPertTimes = dataPertTimes{j, k, :, :};
        currTrueMats = trueMats{j, k, :, :};
        
        save(sprintf('%s/dataLog.mat', currExpPath), 'currDataLog');
        save(sprintf('%s/dataPertLength.mat', currExpPath), 'currDataPertLength');
        save(sprintf('%s/dataPertTimes.mat', currExpPath), 'currDataPertTimes');
        save(sprintf('%s/trueMats.mat', currExpPath), 'currTrueMats');
    end
end

% Save experiment results
save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/tprLog.mat', resultPath), 'tprLog');
save(sprintf('%s/fprLog.mat', resultPath), 'fprLog');
save(sprintf('%s/accLog.mat', resultPath), 'accLog');
save(sprintf('%s/numRerun.mat', resultPath), 'numRerun');


%% Plot Results

forceInd = 1;

% Show number of simulations that were skipped.
figure(1)
imagesc(squeeze(numRerun(:, forceInd, :)))
colorbar
title('Number of Simulations Rerun by Our Analysis')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTickLabel', strengths)
set(gca, 'YTickLabel', networkSizes)
set(gca,'TickLength', [0 0])
set(gca,'YDir','normal')
colormap jet


% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 4);
figure(2)
clims = [0, 1];
imagesc(squeeze(aveAccuracies(:, forceInd, :)), clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
colorbar
title('Average Accuracy over Simulations')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTickLabel', strengths)
set(gca, 'YTickLabel', networkSizes)
set(gca, 'TickLength', [0 0])


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 4);
figure(3)
clims = [0, 1];
imagesc(squeeze(aveTPR(:, forceInd, :)), clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
colorbar
title('Average TPR over Simulations')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTickLabel', strengths)
set(gca, 'YTickLabel', networkSizes)
set(gca, 'TickLength', [0 0])


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 4);
figure(4)
clims = [0, 1];
imagesc(squeeze(aveFPR(:, forceInd, :)), clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
colorbar
title('Average FPR over Simulations')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTickLabel', strengths)
set(gca, 'YTickLabel', networkSizes)
set(gca, 'TickLength', [0 0])
