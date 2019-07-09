clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts')
addpath('../DataScripts/SimulateData')
addpath('../DataScripts/SimulateData/InitFunctions')

expNum = 'VarySizeStrengths';

networkSizes = 20; %2 : 2 : 20;
numSizes = length(networkSizes);

strengths = 1; %1 : 10;
numStrengths = length(strengths);

% Initialize masses, positions, and velocities of oscillators.
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Define time sampling.
deltat = 0.1; % space between time points
endtime = 25;
nobs = round(endtime / deltat); % number of time points (observations)
tSpan = linspace(0, endtime, nobs);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, measParam);

% Specify boundary conditions.
bc = 'fixed';

% Probabilities of network connections.
prob = 0.5;

% Number of matrices to average results over.
numMats = 1; %100;

numTrials = 100;

% Spectral radius threshold for MVGC toolbox.
rhoThresh = 0.995;


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
resultPath = sprintf('%s/GCResults', expPath);
if exist(resultPath, 'dir') == 7
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
    rmdir(resultPath, 's')
end
mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

% Run Granger Causality to infer network connections.
preprocfn = @(data) cos(data);
save(sprintf('%s/expParams.mat', resultPath), 'preprocfn')

predMats = cell(1, numSizes * numStrengths * numMats);
tprLog = nan(1, numSizes * numStrengths * numMats);
fprLog = nan(1, numSizes * numStrengths * numMats);
accLog = nan(1, numSizes * numStrengths * numMats);
numRerun = zeros(1, numSizes * numStrengths * numMats);
diagnosticsLog = nan(numSizes * numStrengths * numMats, 3);

parDataSave = @(fname, noisyData, mat)...
            save(fname, 'noisyData', 'mat');
parResultsSave = @(fname, est, tpr, fpr, acc, diagnostics)...
            save(fname, 'est', 'tpr', 'fpr', 'acc', 'diagnostics');

% Number of parallel processes
M = 12;
c = progress(numSizes * numStrengths * numMats);
for idx = 1 : numSizes * numStrengths * numMats %parfor (idx = 1 : numSizes * numStrengths * numMats, M)
    [j, k, m] = ind2sub([numSizes, numStrengths, numMats], idx);
    fprintf('size: %d, strength: %d\n', j, k)
    
    currExpPath = sprintf('%s/size%d/strength%d/mat%d', expPath, j, k, m);
    if exist(sprintf('%s/dataLog.mat', currExpPath), 'file') ~= 2
        mkdir(currExpPath)
    else
        continue
    end
    
    % Count the number of iterations done by the parfor loop
    c.count();
    
    nvars = networkSizes(j);
    strength = strengths(k);
    
    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);
        
        % Simulate oscillator trajectories.
        forcingFunc = zeros([nvars, nobs]);
        data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
        noisyData = noisefn(data);

        dataObsIdx = true([1, nvars]); % default parameter
        gcSimulationLength = round(max(3, round(22.5 / strength)) / deltat); % number of time points to give to GC
        [est, tableResults] = GrangerBaseExperiment(noisyData(:, 1:gcSimulationLength, :), ...
                mat, preprocfn, dataObsIdx, rhoThresh);
        if isnan(est)
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end
        
        parDataSave(sprintf('%s/dataLog.mat', currExpPath), noisyData, mat);
        parResultsSave(sprintf('%s/results.mat', currExpPath), est, ...
            tableResults.tpr, tableResults.fpr, tableResults.acc, tableResults.diagnostics);
        
        predMats{idx} = est;
        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accLog(idx) = tableResults.acc;
        diagnosticsLog(idx, :) = tableResults.diagnostics;
        break
    end
end

% Reshape data structures
predMats = reshape(predMats, numSizes, numStrengths, numMats);
tprLog = reshape(tprLog, [numSizes, numStrengths, numMats]);
fprLog = reshape(fprLog, [numSizes, numStrengths, numMats]);
accLog = reshape(accLog, [numSizes, numStrengths, numMats]);
diagnosticsLog = reshape(diagnosticsLog, [numSizes, numStrengths, numMats, 3]);
numRerun = sum(reshape(numRerun, [numSizes, numStrengths, numMats]), 3);

% Save experiment results
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', ...
    'accLog', 'diagnosticsLog', 'numRerun');


%% Plot Results

% Show number of simulations that were skipped.
figure(1)
imagesc(reshape(numRerun, [numSizes, numStrengths]))
set(gca,'YDir','normal')
colormap jet
colorbar
title('Number of Simulations Rerun by Our Analysis')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTick', strengths)
set(gca, 'YTick', networkSizes)
%set(gca,'TickLength', [0 0])


% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 3);
figure(2)
clims = [0, 1];
imagesc(reshape(aveAccuracies, [numSizes, numStrengths]), clims)
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
aveTPR = nanmean(tprLog, 3);
figure(3)
clims = [0, 1];
imagesc(reshape(aveTPR, [numSizes, numStrengths]), clims)
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
aveFPR = nanmean(fprLog, 3);
figure(4)
clims = [0, 1];
imagesc(reshape(aveFPR, [numSizes, numStrengths]), clims)
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
