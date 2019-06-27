clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts')
addpath('../DataScripts/SimulateData')
addpath('../DataScripts/SimulateData/InitFunctions')

expNum = 'VaryPertsObs';

% Network size
nvars = 10;

% Connection strengths
numPertsList = 1 : nvars;
numPertsLength = length(numPertsList);

% Forcing magnitudes
numObsList = 1 : nvars;
numObsLength = length(numObsList);

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, measParam);

% Delta t
deltat = 0.01;

% Preprocessing function for data.
preprocfn = @(data) data;

% Amount of time to wait between perturbations
waitTime = 10;

% Probabilities of network connections.
prob = 0.5;

% Coupling strength of network connections.
strength = 50;

% Magnitude of forcing in perturbations.
force = 50;

% Number of matrices to average results over.
numMats = 100;

% Number of experimental trials
numTrials = 1;

% Width of moving variance window
movvarWidth = 10;

% Threshold for mean variance algorithm
meanThresh = 0;

method = 'meanvar';

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
resultPath = sprintf('%s/PertResults', expPath);
if exist(resultPath, 'dir') == 7
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
    rmdir(resultPath, 's')
end
mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

% Run PCI to infer network connections.
predMats = cell(1, numPertsLength * numObsLength * numMats);
tprLog = nan(1, numPertsLength * numObsLength * numMats);
fprLog = nan(1, numPertsLength * numObsLength * numMats);
accLog = nan(1, numPertsLength * numObsLength * numMats);

parsave = @(fname, noisyData, pertIdx, obsIdx, pertLength, pertTimes, mat)...
            save(fname, 'noisyData', 'pertIdx', 'obsIdx', 'pertLength', 'pertTimes', 'mat');

% Number of parallel processes
M = 12;
c = progress(numPertsLength * numObsLength * numMats);
parfor (idx = 1 : numPertsLength * numObsLength * numMats, M)
    [j, k, m] = ind2sub([numPertsLength, numObsLength, numMats], idx);
    fprintf('perts: %d, obs: %d\n', j, k)
    
    % Count the number of iterations done by the parfor loop
    c.count();
    
    numPerts = numPertsList(j);
    numObs = numObsList(k);
    
    if numPerts > numObs
        continue
    end
    
    currExpPath = sprintf('%s/numPerts%d/numObs%d/mat%d', expPath, j, k, m);
    if exist(sprintf('%s/dataLog.mat', currExpPath), 'file') ~= 2
        mkdir(currExpPath)
    else
        continue
    end
    
    % Choose which nodes to observe.
    indsToObserve = randsample(1 : nvars, numObs);
    % Choose which nodes to perturb from the ones you observed.
    pertIdx = randsample(indsToObserve, numPerts);
    
    % Max obsIdx a logical vector.
    obsIdx = false(1, nvars);
    obsIdx(indsToObserve) = true;
    
    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);
        
        endtime = waitTime * (numPerts + 1);
        nobs = round(endtime / deltat);
        tSpan = linspace(0, endtime, nobs);

        % Build up forcing function.
        times = round(linspace(0, nobs, numPerts+2));
        pertTimes = times(2:end-1);
        pertLength = round(waitTime / (4 * deltat));

        forcingFunc = zeros([nvars, nobs]);
        for p=1:numPerts
            forcingFunc(pertIdx(p), pertTimes(p):pertTimes(p)+pertLength) = force;
        end

        % Generate data with forced perturbations.
        data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
        noisyData = noisefn(data);
        
        leftPad = 0;
        rightPad = pertLength;
        truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
        [est, ~, ~, ~, ~, tableResults] = ...
            PerturbationBaseExperiment(noisyData, mat, numTrials, preprocfn, ...
                obsIdx, pertIdx, pertTimes, leftPad, rightPad, method, meanThresh, movvarWidth);
        
        parsave(sprintf('%s/dataLog.mat', currExpPath), noisyData, pertIdx, obsIdx, pertLength, pertTimes, mat);
        
        predMats{idx} = est;
        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accLog(idx) = tableResults.acc;
        break
    end
end

predMats = reshape(predMats, numPertsLength, numObsLength, numMats);
tprLog = reshape(tprLog, [numPertsLength, numObsLength, numMats]);
fprLog = reshape(fprLog, [numPertsLength, numObsLength, numMats]);
accLog = reshape(accLog, [numPertsLength, numObsLength, numMats]);

% Save experiment results
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', 'accLog');


%% Plot Results

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 3);
figure(2)
clims = [0, 1];
imagesc(reshape(aveAccuracies, [numPertsLength, numObsLength]), clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average Accuracy over Simulations')
xlabel('Number of Observations')
ylabel('Number of Perturbations')
set(gca, 'XTick', numObsList)
set(gca, 'YTick', numPertsList)
%set(gca, 'TickLength', [0 0])


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 3);
figure(3)
clims = [0, 1];
imagesc(reshape(aveTPR, [numPertsLength, numObsLength]), clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average TPR over Simulations')
xlabel('Number of Observations')
ylabel('Number of Perturbations')
set(gca, 'XTick', numObsList)
set(gca, 'YTick', numPertsList)
%set(gca, 'TickLength', [0 0])


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 3);
figure(4)
clims = [0, 1];
imagesc(reshape(aveFPR, [numPertsLength, numObsLength]), clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average FPR over Simulations')
xlabel('Number of Observations')
ylabel('Number of Perturbations')
set(gca, 'XTick', numObsList)
set(gca, 'YTick', numPertsList)
%set(gca, 'TickLength', [0 0])
