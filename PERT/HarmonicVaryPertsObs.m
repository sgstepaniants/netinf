clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'PertVarySizeForcingStrengths';

% Network size
nvars = 10;

% Connection strength
strength = 0.1;

% Forcing magnitude
force = 50;

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.3;
cfn = @(n) constfn(n, damping);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, measParam);
preprocfn = @(data) data;

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
