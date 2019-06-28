clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts')
addpath('../DataScripts/SimulateData')
addpath('../DataScripts/SimulateData/InitFunctions')

expNum = 'VaryNoise';

% Network size
nvars = 10;

% Specify noise magnitudes to try
noiseMagnitudes = 0 : 0.2 : 2;
noiseMagnitudesLength = length(noiseMagnitudes);

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

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
numMats = 1;

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
predMats = cell(1, noiseMagnitudesLength * numMats);
tprLog = nan(1, noiseMagnitudesLength * numMats);
fprLog = nan(1, noiseMagnitudesLength * numMats);
accLog = nan(1, noiseMagnitudesLength * numMats);

parsave = @(fname, noisyData, pertIdx, obsIdx, pertLength, pertTimes, mat)...
            save(fname, 'noisyData', 'pertIdx', 'obsIdx', 'pertLength', 'pertTimes', 'mat');

% Number of parallel processes
M = 25;
c = progress(noiseMagnitudesLength * numMats);
parfor (idx = 1 : noiseMagnitudesLength * numMats, M)
    [j, m] = ind2sub([noiseMagnitudesLength, numMats], idx);
    fprintf('noise magnitude: %d\n', j)
    
    currExpPath = sprintf('%s/noise%d/mat%d', expPath, j, m);
    if exist(currExpPath, 'dir') ~= 7
        mkdir(currExpPath)
    end
    
    % Count the number of iterations done by the parfor loop
    c.count();

    noiseMagnitude = noiseMagnitudes(j);

    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);
        
        % Perturb all nodes sequentially.
        pertIdx = 1 : nvars;
        numPerts = length(pertIdx);
        
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
        noisyData = WhiteGaussianNoise(data, noiseMagnitude);
        
        obsIdx = true([1, nvars]);
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


predMats = reshape(predMats, noiseMagnitudesLength, numMats);
tprLog = reshape(tprLog, [noiseMagnitudesLength, numMats]);
fprLog = reshape(fprLog, [noiseMagnitudesLength, numMats]);
accLog = reshape(accLog, [noiseMagnitudesLength, numMats]);

% Save experiment results
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', 'accLog');


%% Plot Results

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 2);
figure(1)
plot(noiseMagnitudes, aveAccuracies)
xlabel('Magnitude of Noise')
ylabel('Average Accuracy over Simulations')


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 4);
figure(2)
plot(noiseMagnitudes, aveTPR)
xlabel('Magnitude of Noise')
ylabel('Average TPR over Simulations')


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 4);
figure(3)
plot(noiseMagnitudes, aveFPR)
xlabel('Magnitude of Noise')
ylabel('Average FPR over Simulations')
