clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

% Network size
nvars = 2;

% Experiment name
expNum = sprintf('VaryStrengthsProbs_Size%d', nvars);

% Initialize positions and frequencies of oscillators.
pfn = @(n) randfn(n, -0.5, 0.5);
wfn = @(n) randfn(n, -1, 1);

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

% Probabilities of network connections to try.
probs = 0.1:0.1:1;
numProbs = length(probs);

% Connection strengths to try.
strengths = 1:1:10;
numStrengths = length(strengths);

% Number of matrices to try for each probability and connection strength
% combination.
numMats = 1;

% Number of simulation trials.
numTrials = 10;

% Spectral radius threshold for MVGC toolbox.
rhoThresh = 0.995;


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
dataLog = nan(nvars, nobs, numTrials, numProbs * numStrengths * numMats);
trueMats = nan(nvars, nvars, numProbs * numStrengths * numMats);

% Run Granger Causality to infer network connections.
freq = 1;
preprocfn = @(data) cos(data);
save(sprintf('%s/expParams.mat', resultPath), 'freq', 'preprocfn')

predMats = nan(nvars, nvars, numProbs * numStrengths * numMats);
numRerun = zeros(1, numProbs * numStrengths * numMats);
tprLog = nan(1, numProbs * numStrengths * numMats);
fprLog = nan(1, numProbs * numStrengths * numMats);
accuracyLog = nan(1, numProbs * numStrengths * numMats);
diagnosticsLog = nan(numProbs * numStrengths * numMats, 3);

% Number of parallel processes
M = 12;
c = progress(numProbs * numStrengths * numMats);
parfor (idx = 1 : numProbs * numStrengths * numMats, M)
    [j, k, l] = ind2sub([numProbs, numStrengths, numMats], idx);
    prob = probs(j);
    strength = strengths(k);
    
    c.count();
    
    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);
        
        % Specify forcing function for oscillators.
        forcingFunc = zeros([nvars, nobs]);

        % Simulate oscillator trajectories.
        data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
        noisyData = noisefn(data);

        dataObsIdx = true([1, nvars]); % default parameter
        [est, tableResults] = GrangerBaseExperiment(noisyData, ...
                mat, preprocfn, freq, '', dataObsIdx, rhoThresh);
        if isnan(est)
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end

        dataLog(:, :, :, idx) = noisyData;
        trueMats(:, :, idx) = mat;
        predMats(:, :, idx) = est;

        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accuracyLog(idx) = tableResults.acc;
        diagnosticsLog(idx, :) = tableResults.diagnostics;
        break
    end
end

% Reshape data structures
dataLog = reshape(dataLog, [nvars, nobs, numTrials, numProbs, numStrengths, numMats]);
trueMats = reshape(trueMats, [nvars, nvars, numProbs, numStrengths, numMats]);
predMats = reshape(predMats, [nvars, nvars, numProbs, numStrengths, numMats]);
numRerun = sum(reshape(numRerun, [numProbs, numStrengths, numMats]), 3);
tprLog = reshape(tprLog, [numProbs, numStrengths, numMats]);
fprLog = reshape(fprLog, [numProbs, numStrengths, numMats]);
accuracyLog = reshape(accuracyLog, [numProbs, numStrengths, numMats]);
diagnosticsLog = reshape(diagnosticsLog, [numProbs, numStrengths, numMats, 3]);

% Save experiment simulated data and connectivity matrices.
save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/trueMats.mat', expPath), 'trueMats');

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
imagesc(numRerun)
colorbar
title('Number of Simulations Rerun by Our Analysis')
xlabel('Connection Strength')
ylabel('Connection Probability')
set(gca, 'XTickLabel', strengths)
set(gca, 'YTickLabel', probs)
set(gca,'TickLength',[0 0])
set(gca,'YDir','normal')
colormap jet


% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accuracyLog, 3);
figure(2)
clims = [0, 1];
imagesc(aveAccuracies, clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
colorbar
title('Average Accuracy over Simulations')
xlabel('Connection Strength')
ylabel('Connection Probability')
set(gca, 'XTickLabel', strengths)
set(gca, 'YTickLabel', probs)
set(gca,'TickLength',[0 0])


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 3);
figure(3)
clims = [0, 1];
imagesc(aveTPR, clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
colorbar
title('Average TPR over Simulations')
xlabel('Connection Strength')
ylabel('Connection Probability')
set(gca, 'XTickLabel', strengths)
set(gca, 'YTickLabel', probs)
set(gca,'TickLength',[0 0])


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 3);
figure(4)
clims = [0, 1];
imagesc(aveFPR, clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
colorbar
title('Average FPR over Simulations')
xlabel('Connection Strength')
ylabel('Connection Probability')
set(gca, 'XTickLabel', strengths)
set(gca, 'YTickLabel', probs)
set(gca,'TickLength',[0 0])