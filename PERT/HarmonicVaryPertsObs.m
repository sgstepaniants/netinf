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

% Delta t
deltat = 0.1;

% Specify boundary conditions.
bc = 'fixed';

% Preprocessing function for data.
preprocfn = @(data) data;

% Probabilities of network connections.
prob = 0.5;

% Coupling strength of network connections.
strength = 0.1;

% Magnitude of forcing in perturbations.
force = 50;

% Number of matrices to average results over.
numMats = 1;

% Number of experimental trials
numTrials = 1;

method = 'corr';

% Threshold for correlation algorithm
corrThresh = 0.2;

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../HarmonicExperiments/%s', expName);
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
numRerun = zeros(1, numPertsLength * numObsLength * numMats);

% Number of parallel processes
M = 12;
c = progress(numPertsLength * numObsLength * numMats);
for idx = 1 : numPertsLength * numObsLength * numMats %parfor (idx = 1 : numPertsLength * numObsLength * numMats, M)
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
        K = MakeNetworkTriDiag(nvars+2, false);
        K(2:nvars+1, 2:nvars+1) = mat;
        K = strength * K;

        % If this adjacency matrix is bad, make a new simulation.
        [disconnectedNodes, amplitudes, waitTime] = checkHarmonicMat(K, damping, force);
        if waitTime > 500 || ~isempty(disconnectedNodes) || any(amplitudes > -0.00001)
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
        for p=1:numPerts
            forcingFunc(pertIdx(p), pertTimes(p):pertTimes(p)+pertLength) = force;
        end

        % Generate data with forced perturbations.
        data = GenerateHarmonicData(nvars, tSpan, numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
        noisyData = noisefn(data);
        
        leftPad = 100;
        rightPad = pertLength;
        [est, ~, ~, ~, ~, tableResults] = ...
            PerturbationBaseExperiment(noisyData, mat, numTrials, preprocfn, ...
                obsIdx, pertIdx, pertTimes, leftPad, rightPad, method, corrThresh);
        
        parSave.parPertDataSave(sprintf('%s/dataLog.mat', currExpPath), noisyData, pertIdx, obsIdx, pertLength, pertTimes, mat, K);
        
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
numRerun = sum(reshape(numRerun, [numPertsLength, numObsLength, numMats]), 3);

% Save experiment results
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', 'accLog', 'numRerun');


%% Plot Results

% Show number of simulations that were skipped.
figure(1)
imagesc(reshape(numRerun, [numPertsLength, numObsLength]))
set(gca,'YDir','normal')
colormap jet
colorbar
title('Number over Simulations Rerun by Analysis')
xlabel('Number of Observations')
ylabel('Number of Perturbations')
set(gca, 'XTick', numObsList)
set(gca, 'YTick', numPertsList)
%set(gca, 'TickLength', [0 0])


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
