clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')
addpath('../kmeans_opt/')

nvars = 5;
prob = 0.5;
spring = 0.1;
damping = 0.3;
pertForce = 30;

expNum = sprintf('EXPVaryPertsObs(nvars%d_prob%.2f_spring%.2f_damping%.2f_pertf%.2f)', nvars, prob, spring, damping, pertForce);
dataPath = sprintf('../HarmonicExperiments/%s', expNum);

if exist(dataPath, 'dir') ~= 7
    error('Data not found: %s', dataPath)
end

outputPath = sprintf('%s/PertResults', dataPath);
if exist(outputPath, 'dir') ~= 7
    mkdir(outputPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', outputPath),'s');
    if upper(m) == 'N'
       return
    end
end

% Threshold for correlation algorithm
corrThresh = 0.5;

% Padding for window
pad = 1;

% Method for building probability matrix
method = 'corr';

% How often to save experiment results
freq = 1;

% Function to preprocess data
preprocfn = @(data) data;

%% Evaluate Algorithm on Data for Varying Numbers of Perturbations and Observations
load(sprintf('%s/params.mat', dataPath))

predMats = nan(nvars, nvars, numTrials, numMats, nvars, nvars);
tprLog = nan(nvars, nvars, numTrials, numMats);
fprLog = nan(nvars, nvars, numTrials, numMats);
accuracyLog = nan(nvars, nvars, numTrials, numMats);

fprintf('Run Algorithm:\n')
for numObs = nvars:-1:1
    numObs
    for numPerts = nvars:-1:1
        numPerts
        
        try
            load(sprintf('%s/numobs%d/numperts%d/dataLog.mat', dataPath, numObs, numPerts), 'dataLog');
            load(sprintf('%s/numobs%d/numperts%d/trueMats.mat', dataPath, numObs, numPerts), 'trueMats');
            load(sprintf('%s/numobs%d/numperts%d/dataObsIdx.mat', dataPath, numObs, numPerts), 'dataObsIdx');
            load(sprintf('%s/numobs%d/numperts%d/dataPertIdx.mat', dataPath, numObs, numPerts), 'dataPertIdx');
            load(sprintf('%s/numobs%d/numperts%d/dataPertTimes.mat', dataPath, numObs, numPerts), 'dataPertTimes');
            load(sprintf('%s/numobs%d/numperts%d/dataPertLength.mat', dataPath, numObs, numPerts), 'dataPertLength');
        catch
            continue
        end
        
        % Run Perturbation Inference to infer network connections.
        [predMats(:, :, :, :, numObs, numPerts), tableResults] = ...
            PerturbationBaseExperiment(dataLog, trueMats, numTrials, preprocfn, ...
                dataObsIdx, dataPertIdx, dataPertTimes, dataPertLength, ...
                method, corrThresh, pad, 0, freq, outputPath);
        
        tprLog(numObs, numPerts, :, :) = tableResults.tpr;
        fprLog(numObs, numPerts, :, :) = tableResults.fpr;
        accuracyLog(numObs, numPerts, :, :) = tableResults.acc;
        
        save(sprintf('%s/predMats.mat', outputPath), 'predMats');
        save(sprintf('%s/tprLog.mat', outputPath), 'tprLog');
        save(sprintf('%s/fprLog.mat', outputPath), 'fprLog');
        save(sprintf('%s/accuracyLog.mat', outputPath), 'accuracyLog');
    end
end

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(nanmean(accuracyLog, 4), 3);
figure(1)
clims = [0, 1];
imagesc(aveAccuracies, clims)
colorbar
title('Average Accuracy over Simulations')
xlabel('Number of Perturbed Nodes')
ylabel('Number of Observed Nodes')

% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(nanmean(tprLog, 4), 3);
figure(2)
clims = [0, 1];
imagesc(aveTPR, clims)
colorbar
title('Average TPR over Simulations')
xlabel('Number of Perturbed Nodes')
ylabel('Number of Observed Nodes')

% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(nanmean(fprLog, 4), 3);
figure(3)
clims = [0, 1];
imagesc(aveFPR, clims)
colorbar
title('Average FPR over Simulations')
xlabel('Number of Perturbed Nodes')
ylabel('Number of Observed Nodes')
