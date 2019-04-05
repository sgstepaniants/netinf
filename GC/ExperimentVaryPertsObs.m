clear all; close all; clc;
addpath('../DataScripts/SimulateData/')

nvars = 5;
prob = 0.5;
spring = 0.1;
damping = 0.3;
pertForce = 30;

expNum = sprintf('EXPVaryPertsObs(nvars%d_prob%.2f_spring%.2f_damping%.2f_pertf%.2f)', nvars, prob, spring, damping, pertForce);
expPath = sprintf('../HarmonicExperiments/%s', expNum);

if exist(expPath, 'dir') ~= 7
    error('Data not found: %s', expPath)
end

resultPath = sprintf('%s/GCResults', expPath);
if exist(resultPath, 'dir') ~= 7
    mkdir(resultPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
end

% Number of repetitions of GC
reps = 1;

% How often to save GC results
freq = 1;

% Function to preprocess data
preprocfn = @(data) data;
%preprocfn = @(data) mydetrend(data);


%% Evaluate Algorithm on Data for Varying Numbers of Perturbations and Observations
load(sprintf('%s/params.mat', expPath))

predMats = nan(nvars, nvars, numMats, nvars, nvars);
normLog = nan(nvars, nvars, numMats);
tprLog = nan(nvars, nvars, numMats);
fprLog = nan(nvars, nvars, numMats);
accuracyLog = nan(nvars, nvars, numMats);
tableResultsLog = repmat(struct([]), [nvars, nvars]);

fprintf('Run Algorithm:\n')
for numObs = nvars:-1:2
    numObs
    for numPerts = nvars:-1:1
        numPerts
        
        try
            load(sprintf('%s/numobs%d/numperts%d/dataLog.mat', expPath, numObs, numPerts), 'dataLog');
            load(sprintf('%s/numobs%d/numperts%d/trueMats.mat', expPath, numObs, numPerts), 'trueMats');
            load(sprintf('%s/numobs%d/numperts%d/dataObsIdx.mat', expPath, numObs, numPerts), 'dataObsIdx');
            %load(sprintf('%s/numobs%d/numperts%d/dataPertIdx.mat', expPath, numObs, numPerts), 'dataPertIdx');
            %load(sprintf('%s/numobs%d/numperts%d/dataPertTimes.mat', expPath, numObs, numPerts), 'dataPertTimes');
            %load(sprintf('%s/numobs%d/numperts%d/dataPertLength.mat', expPath, numObs, numPerts), 'dataPertLength');
        catch
            continue
        end
        
        % Run Granger Causality to infer network connections.
        [predMats(:, :, :, numObs, numPerts), ~, tableResults] = ...
            GrangerBaseExperiment(dataLog, trueMats, numTrials, reps, preprocfn, freq, resultPath, dataObsIdx);
        
        tableResultsLog(numObs, numPerts).norm = tableResults.norm;
        tableResultsLog(numObs, numPerts).normVoting = tableResults.normVoting;
        tableResultsLog(numObs, numPerts).tpr = tableResults.tpr;
        tableResultsLog(numObs, numPerts).tprVoting = tableResults.tprVoting;
        tableResultsLog(numObs, numPerts).fpr = tableResults.fpr;
        tableResultsLog(numObs, numPerts).fprVoting = tableResults.fprVoting;
        tableResultsLog(numObs, numPerts).acc = tableResults.acc;
        tableResultsLog(numObs, numPerts).accVoting = tableResults.accVoting;
        tableResultsLog(numObs, numPerts).diagnostics = tableResults.diagnostics;
        
        normLog(numObs, numPerts, :) = tableResults.normVoting;
        tprLog(numObs, numPerts, :) = tableResults.tprVoting;
        fprLog(numObs, numPerts, :) = tableResults.fprVoting;
        accuracyLog(numObs, numPerts, :) = tableResults.accVoting;
        
        save(sprintf('%s/predMats.mat', resultPath), 'predMats');
        save(sprintf('%s/tableResultsLog.mat', resultPath), 'tableResultsLog');
        save(sprintf('%s/normLog.mat', resultPath), 'normLog');
        save(sprintf('%s/tprLog.mat', resultPath), 'tprLog');
        save(sprintf('%s/fprLog.mat', resultPath), 'fprLog');
        save(sprintf('%s/accuracyLog.mat', resultPath), 'accuracyLog');
    end
end

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accuracyLog, 3);
figure(1)
clims = [0, 1];
imagesc(aveAccuracies, clims)
colorbar
title('Average Accuracy over Simulations')
xlabel('Number of Perturbed Nodes')
ylabel('Number of Observed Nodes')

% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 3);
figure(2)
clims = [0, 1];
imagesc(aveTPR, clims)
colorbar
title('Average TPR over Simulations')
xlabel('Number of Perturbed Nodes')
ylabel('Number of Observed Nodes')

% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 3);
figure(3)
clims = [0, 1];
imagesc(aveFPR, clims)
colorbar
title('Average FPR over Simulations')
xlabel('Number of Perturbed Nodes')
ylabel('Number of Observed Nodes')
