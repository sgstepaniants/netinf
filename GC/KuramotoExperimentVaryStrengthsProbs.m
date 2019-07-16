clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

% Network size
nvars = 5;

% Experiment name
expNum = sprintf('VaryStrengthsProbs_Size%d', nvars);

% Probabilities of network connections to try.
probs = 0.5; %0 : 0.1 : 1;
numProbs = length(probs);

% Connection strengths to try.
strengths = 5; %1 : 1 : 10;
numStrengths = length(strengths);

% Initialize positions and frequencies of oscillators.
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Define time sampling.
deltat = 0.1; % space between time points

endtime = 25;
nobs = round(endtime / deltat); % number of time points (observations)
tSpan = linspace(0, endtime, nobs);

% Specify forcing function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

% Number of matrices to try for each probability and connection strength
% combination.
numMats = 1; %100;

% Number of simulation trials.
numTrials = 100;

% Spectral radius threshold for MVGC toolbox.
rhoThresh = 0.995;


% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);
%if exist(expPath, 'dir') == 7
%    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
%    if upper(m) == 'N'
%        return
%    end
%    rmdir(expPath, 's')
%end
%mkdir(expPath)

% Save experiment parameters.
%save(sprintf('%s/params.mat', expPath));

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/GCResults', expPath);
%if exist(resultPath, 'dir') == 7
%    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
%    if upper(m) == 'N'
%       return
%    end
%    rmdir(resultPath, 's')
%end
%mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

% Run Granger Causality to infer network connections.
preprocfn = @(data) cos(data);
save(sprintf('%s/expParams.mat', resultPath), 'preprocfn')

predMats = nan(nvars, nvars, numProbs * numStrengths * numMats);
tprLog = nan(1, numProbs * numStrengths * numMats);
fprLog = nan(1, numProbs * numStrengths * numMats);
accLog = nan(1, numProbs * numStrengths * numMats);
numRerun = zeros(1, numProbs * numStrengths * numMats);
diagnosticsLog = nan(numProbs * numStrengths * numMats, 3);

parDataSave = @(fname, noisyData, mat)...
            save(fname, 'noisyData', 'mat');
parResultsSave = @(fname, est, tpr, fpr, acc, diagnostics)...
            save(fname, 'est', 'tpr', 'fpr', 'acc', 'diagnostics');

% Number of parallel processes
M = 12;
c = progress(numProbs * numStrengths * numMats);
for idx = 1 : numProbs * numStrengths * numMats %parfor (idx = 1 : numProbs * numStrengths * numMats, M)
    [j, k, m] = ind2sub([numProbs, numStrengths, numMats], idx);
    fprintf('prob: %d, strength: %d\n', j, k)
    
    c.count();
    
    currExpPath = sprintf('%s/prob%d/strength%d/mat%d', expPath, j, k, m);
    if exist(sprintf('%s/dataLog.mat', currExpPath), 'file') ~= 2
        mkdir(currExpPath)
    end
    
    prob = probs(j);
    strength = strengths(k);
    
    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);
        
        % Simulate oscillator trajectories.
        data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
        noisyData = noisefn(data);

        dataObsIdx = true([1, nvars]);
        % number of time points to give to GC
        gcSimulationLength = round(max(3, min(round(4.5 * nvars / strength), 25)) / deltat);
        [est, tableResults] = GrangerBaseExperiment(noisyData(:, 1:gcSimulationLength, :), ...
                mat, preprocfn, dataObsIdx, rhoThresh);
        if isnan(est)
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end
        
        %parDataSave(sprintf('%s/dataLog.mat', currExpPath), noisyData, mat);
        %parResultsSave(sprintf('%s/results.mat', currExpPath), est, ...
        %    tableResults.tpr, tableResults.fpr, tableResults.acc, tableResults.diagnostics);
        
        predMats(:, :, idx) = est;
        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accLog(idx) = tableResults.acc;
        diagnosticsLog(idx, :) = tableResults.diagnostics;
        break
    end
end

% Reshape data structures
predMats = reshape(predMats, [nvars, nvars, numProbs, numStrengths, numMats]);
tprLog = reshape(tprLog, [numProbs, numStrengths, numMats]);
fprLog = reshape(fprLog, [numProbs, numStrengths, numMats]);
accLog = reshape(accLog, [numProbs, numStrengths, numMats]);
diagnosticsLog = reshape(diagnosticsLog, [numProbs, numStrengths, numMats, 3]);
numRerun = sum(reshape(numRerun, [numProbs, numStrengths, numMats]), 3);

%save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', ...
%    'accLog', 'diagnosticsLog', 'numRerun');


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
set(gca,'YDir','normal')
colormap jet


% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 3);
figure(2)
clims = [0, 1];
imagesc(aveAccuracies, clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
%colorbar
%title('Average Accuracy over Simulations')
%xlabel('Connection Strength')
%ylabel('Connection Probability')
%set(gca, 'XTickLabel', strengths)
%set(gca, 'YTickLabel', probs)


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
%colorbar
%title('Average TPR over Simulations')
%xlabel('Connection Strength')
%ylabel('Connection Probability')
%set(gca, 'XTickLabel', strengths)
%set(gca, 'YTickLabel', probs)


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
%colorbar
%title('Average FPR over Simulations')
%xlabel('Connection Strength')
%ylabel('Connection Probability')
%set(gca, 'XTickLabel', strengths)
%set(gca, 'YTickLabel', probs)
