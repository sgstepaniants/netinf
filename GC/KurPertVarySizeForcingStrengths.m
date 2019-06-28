clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')

expNum = 'PertVarySizeForcingStrengths';

% Preprocessing function for data.
preprocfn = @(data) downsample(cos(data).', 20).';

% Spectral radius threshold for MVGC toolbox.
rhoThresh = 1;

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);

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
load(sprintf('%s/params.mat', expPath), '-regexp', '^(?!expNum$|expName$|expPath$|resultPath$|preprocfn$).')

% Run GC to infer network connections.
predMats = cell(1, numSizes * numForces * numStrengths * numMats);
tprLog = nan(1, numSizes * numForces * numStrengths * numMats);
fprLog = nan(1, numSizes * numForces * numStrengths * numMats);
accLog = nan(1, numSizes * numForces * numStrengths * numMats);
diagnosticsLog = nan(numSizes * numForces * numStrengths * numMats, 3);
numRerun = zeros(1, numSizes * numForces * numStrengths * numMats);

% Number of parallel processes
M = 25;
c = progress(numSizes * numForces * numStrengths * numMats);
for idx = 1 : numSizes * numForces * numStrengths * numMats %parfor (idx = 1 : numSizes * numForces * numStrengths * numMats, M)
    [j, k, l, m] = ind2sub([numSizes, numForces, numStrengths, numMats], idx);
    fprintf('size: %d, force: %d, strength: %d\n', j, k, l)
    
    currExpPath = sprintf('%s/size%d/force%d/strength%d/mat%d', expPath, j, k, l, m);
    if exist(currExpPath, 'dir') ~= 7
        mkdir(currExpPath)
    end
    
    % Count the number of iterations done by the parfor loop
    c.count();

    while true
        dataLog = load(sprintf('%s/dataLog.mat', currExpPath));
        
        data = dataLog.noisyData;
        mat = dataLog.mat;
        nvars = size(mat, 1);
        dataObsIdx = true([1, nvars]); % default parameter
        [est, tableResults] = GrangerBaseExperiment(data, ...
                mat, preprocfn, dataObsIdx, rhoThresh);
        if isnan(est)
            numRerun(idx) = numRerun(idx) + 1;
        end
        
        predMats{idx} = est;
        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accLog(idx) = tableResults.acc;
        diagnosticsLog(idx, :) = tableResults.diagnostics;
        break
    end
end

% Reshape data structures
predMats = reshape(predMats, numSizes, numForces, numStrengths, numMats);
tprLog = reshape(tprLog, [numSizes, numForces, numStrengths, numMats]);
fprLog = reshape(fprLog, [numSizes, numForces, numStrengths, numMats]);
accLog = reshape(accLog, [numSizes, numForces, numStrengths, numMats]);
diagnosticsLog = reshape(diagnosticsLog, [numSizes, numForces, numStrengths, numMats, 3]);
numRerun = sum(reshape(numRerun, [numSizes, numForces, numStrengths, numMats]), 4);

% Save experiment results
save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/tprLog.mat', resultPath), 'tprLog');
save(sprintf('%s/fprLog.mat', resultPath), 'fprLog');
save(sprintf('%s/accLog.mat', resultPath), 'accLog');
save(sprintf('%s/diagnosticsLog.mat', resultPath), 'diagnosticsLog');
save(sprintf('%s/numRerun.mat', resultPath), 'numRerun');


%% Plot Results

forceInd = 1;

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
