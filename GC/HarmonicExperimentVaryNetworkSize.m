clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'PaperVaryNetworkSize';

networkSizes = 2:4;
numSizes = length(networkSizes);

% Initialize masses, positions, and velocities of oscillators.
mfn = @(n) constfn(n, 1); %randfn(n, 0.7, 1.3);
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) zeros([n, 1]);

% Specify the damping constant.
damping = 0.2;
cfn = @(n) constfn(n, damping);

% Define time sampling.
deltat = 0.1; % space between time points
endtime = 25;
nobs = round(endtime / deltat); % number of time points (observations)
tSpan = linspace(0, endtime, nobs);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

% Specify boundary conditions.
bc = 'fixed';

% Probabilities of network connections.
prob = 0.5;

% Connection strength.
strength = 1;

% Number of matrices to average results over.
numMats = 3;

% Number of simulation trials per repetition.
numTrials = 10;

% Spectral radius threshold for MVGC toolbox.
rhoThresh = 0.995;


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
Ks = cell(1, numSizes * numMats);

% Run Granger Causality to infer network connections.
freq = 1;
preprocfn = @(data) standardize(data);
save(sprintf('%s/expParams.mat', resultPath), 'freq', 'preprocfn')

predMats = cell(1, numSizes * numMats);
numRerun = zeros(1, numSizes * numMats);
tprLog = nan(1, numSizes * numMats);
fprLog = nan(1, numSizes * numMats);
accuracyLog = nan(1, numSizes * numMats);
diagnosticsLog = nan(numSizes * numMats, 3);

% Number of parallel processes
M = 12;
parfor (idx = 1 : numSizes * numMats, M)
    [j, l] = ind2sub([numSizes, numMats], idx);
    fprintf('size: %d, mat: %d\n', j, l)
    
    nvars = networkSizes(j);
    
    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);
        K = MakeNetworkTriDiag(nvars+2, false);
        K(2:nvars+1, 2:nvars+1) = mat;
        K = strength * K;
        
        % If any nodes in the network are not connected to the walls or
        % the eigenvalues of the system have positive real parts, don't
        % use this network.
        [~, amplitudes] = checkHarmonicMat(K, damping);
        if any(amplitudes > 0)
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end
        
        % Specify forcing function for oscillators.
        forcingFunc = zeros([nvars, nobs]);

        % Simulate oscillator trajectories.
        fprintf('Generate Data\n')
        data = GenerateHarmonicData(nvars, tSpan, ...
                numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
        noisyData = noisefn(data);

        dataObsIdx = true([1, nvars]); % default parameter
        fprintf('Run Experiment\n')
        [est, tableResults] = GrangerBaseExperiment(noisyData, ...
                mat, preprocfn, freq, '', dataObsIdx, rhoThresh);
        if isnan(est)
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end
        
        dataLog{idx} = noisyData;
        trueMats{idx} = mat;
        Ks{idx} = K;
        predMats{idx} = est;

        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accuracyLog(idx) = tableResults.acc;
        diagnosticsLog(idx, :) = tableResults.diagnostics;

        l = l + 1;
        break
    end
end

% Reshape data structures
dataLog = reshape(dataLog, numSizes, numMats);
trueMats = reshape(trueMats, numSizes, numMats);
Ks = reshape(Ks, numSizes, numMats);
predMats = reshape(predMats, numSizes, numMats);
numRerun = sum(reshape(numRerun, [numSizes, numMats]), 2);
tprLog = reshape(tprLog, [numSizes, numMats]);
fprLog = reshape(fprLog, [numSizes, numMats]);
accuracyLog = reshape(accuracyLog, [numSizes, numMats]);
diagnosticsLog = reshape(diagnosticsLog, [numSizes, numMats, 3]);

% Save experiment simulated data and connectivity matrices.
save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/trueMats.mat', expPath), 'trueMats');
save(sprintf('%s/Ks.mat', expPath), 'Ks');

% Save experiment results
save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/tprLog.mat', resultPath), 'tprLog');
save(sprintf('%s/fprLog.mat', resultPath), 'fprLog');
save(sprintf('%s/accuracyLog.mat', resultPath), 'accuracyLog');
save(sprintf('%s/diagnosticsLog.mat', resultPath), 'diagnosticsLog');
save(sprintf('%s/numRerun.mat', resultPath), 'numRerun');


%% Plot Results

% Show number of simulations that were skipped.
figure(1)
plot(networkSizes, numRerun)
title('Number of Simulations Rerun by Our Analysis')
xlabel('Network Size')
ylabel('Number of Times Experiment Rerun')

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accuracyLog, 2);
figure(2)
plot(networkSizes, aveAccuracies)
title('Average Accuracy over Simulations')
xlabel('Network Size')
ylabel('Accuracy')

% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 2);
figure(3)
plot(networkSizes, aveTPR)
title('Average TPR over Simulations')
xlabel('Network Size')
ylabel('TPR')

% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 2);
figure(4)
plot(networkSizes, aveFPR)
title('Average FPR over Simulations')
xlabel('Network Size')
ylabel('FPR')
