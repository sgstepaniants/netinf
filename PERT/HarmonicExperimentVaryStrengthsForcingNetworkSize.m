clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'PertVaryStrengthsNetworkSize';

% Network sizes
networkSizes = 2 : 20;
numSizes = length(networkSizes);

% Connection strengths
strengths = 0.1 : 0.1 : 1;
numStrengths = length(strengths);

% Forcing magnitudes
forces = 10 : 10 : 50;
numForces = length(forces);

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.3;
cfn = @(n) constfn(n, damping);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

% Delta t
deltat = 0.1;

% Specify boundary conditions.
bc = 'fixed';

% Probabilities of network connections.
prob = 0.5;

% Number of matrices to average results over.
numMats = 100;

% Number of experimental trials
numTrials = 1;

method = 'corr';

% Threshold for correlation algorithm
corrThresh = 0.2;

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../HarmonicExperiments/%s', expName);
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

%% Generate Data and Run Granger Causality Experiments

% Create random connectivity matrices and simulate oscillator trajectories.
dataLog = cell(1, numSizes * numStrengths * numMats);
trueMats = cell(1, numSizes * numStrengths * numMats);
Ks = cell(1, numSizes * numStrengths * numForces * numMats);

% Run PCI to infer network connections.
preprocfn = @(data) data - repmat(uniffn(size(data, 1), -0.5, 0.5, bc), [1, size(data, 2)]);
save(sprintf('%s/expParams.mat', resultPath), 'preprocfn')

numRerun = zeros(1, numSizes * numStrengths * numForces * numMats);
predMats = cell(1, numSizes * numStrengths * numForces * numMats);
tprLog = nan(1, numSizes * numStrengths * numForces * numMats);
fprLog = nan(1, numSizes * numStrengths * numForces * numMats);
accuracyLog = nan(1, numSizes * numStrengths * numForces * numMats);

% Number of parallel processes
M = 25;
c = progress(numSizes * numStrengths * numForces * numMats);
parfor (idx = 1 : numSizes * numStrengths * numForces * numMats, M)
    [j, k, l, m] = ind2sub([numSizes, numStrengths, numForces, numMats], idx);
    fprintf('size: %d, strength: %d, mat: %d\n', j, k, l)

    % Count the number of iterations done by the parfor loop
    c.count();

    nvars = networkSizes(j);
    strength = strengths(k);
    force = forces(l);

    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);
        K = MakeNetworkTriDiag(nvars+2, false);
        K(2:nvars+1, 2:nvars+1) = mat;
        K = strength * K;

        % Perturb all nodes sequentially.
        pertIdx = 1 : nvars;
        numPerts = length(pertIdx);

        % If this adjacency matrix is bad, make a new simulation.
        [disconnectedNodes, amplitudes, waitTime] = checkHarmonicMat(K, damping, force);
        if ~isempty(disconnectedNodes) || any(amplitudes > -0.00001) || waitTime > 500
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end

        endtime = waitTime * (numPerts + 1);
        nobs = round(endtime / deltat);
        tSpan = linspace(0, endtime, nobs);

        % Build up forcing function.
        times = round(linspace(0, nobs, numPerts+2));
        pertTimes = times(2:end-1);
        pertLength = round(waitTime / (4 * deltat));

        forcingFunc = zeros([nvars, nobs]);
        for k=1:numPerts
            forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = force;
        end

        % Generate data with forced perturbations.
        data = GenerateHarmonicData(nvars, tSpan, numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
        noisyData = noisefn(data);

        obsIdx = true([1, nvars]); % default parameter
        leftPad = 100;
        rightPad = pertLength;
        [est, predMatsHist, AprobHist, truePertOrders, predPertOrders, tableResults] = ...
            PerturbationBaseExperiment(noisyData, mat, numTrials, preprocfn, ...
                obsIdx, pertIdx, pertTimes, leftPad, rightPad, method, corrThresh);

        dataLog{idx} = noisyData;
        trueMats{idx} = mat;
        Ks{idx} = K;
        predMats{idx} = est;

        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accuracyLog(idx) = tableResults.acc;
        break
    end
end

% Reshape data structures
dataLog = reshape(dataLog, numSizes, numStrengths, numForces, numMats);
numRerun = sum(reshape(numRerun, [numSizes, numStrengths, numForces, numMats]), 4);
trueMats = reshape(trueMats, numSizes, numStrengths, numForces, numMats);
Ks = reshape(Ks, numSizes, numStrengths, numForces, numMats);
predMats = reshape(predMats, numSizes, numStrengths, numForces, numMats);
tprLog = reshape(tprLog, [numSizes, numStrengths, numForces, numMats]);
fprLog = reshape(fprLog, [numSizes, numStrengths, numForces, numMats]);
accuracyLog = reshape(accuracyLog, [numSizes, numStrengths, numForces, numMats]);

% Save experiment simulated data and connectivity matrices.
save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/trueMats.mat', expPath), 'trueMats');
save(sprintf('%s/Ks.mat', expPath), 'Ks');

% Save experiment results
save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/tprLog.mat', resultPath), 'tprLog');
save(sprintf('%s/fprLog.mat', resultPath), 'fprLog');
save(sprintf('%s/accuracyLog.mat', resultPath), 'accuracyLog');
save(sprintf('%s/numRerun.mat', resultPath), 'numRerun');


%% Plot Results

forceInd = 1;

% Show number of simulations that were skipped.
figure(1)
imagesc(numRerun(:, :, forceInd))
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
aveAccuracies = nanmean(accuracyLog, 4);
figure(2)
clims = [0, 1];
imagesc(aveAccuracies(:, :, forceInd), clims)
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
imagesc(aveTPR(:, :, forceInd), clims)
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
imagesc(aveFPR(:, :, forceInd), clims)
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
