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

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Specify the damping constant.
damping = 0.3;
cfn = @(n) constfn(n, damping);

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
numMats = 100;

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
predMats = cell(1, noiseMagnitudesLength * numMats);
tprLog = nan(1, noiseMagnitudesLength * numMats);
fprLog = nan(1, noiseMagnitudesLength * numMats);
accLog = nan(1, noiseMagnitudesLength * numMats);
numRerun = zeros(1, noiseMagnitudesLength * numMats);

parsave = @(fname, noisyData, pertLength, pertTimes, mat, K)...
            save(fname, 'noisyData', 'pertLength', 'pertTimes', 'mat', 'K');

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
        K = MakeNetworkTriDiag(nvars+2, false);
        K(2:nvars+1, 2:nvars+1) = mat;
        K = strength * K;

        % Perturb all nodes sequentially.
        pertIdx = 1 : nvars;
        numPerts = length(pertIdx);

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
        noisyData = WhiteGaussianNoise(data, noiseMagnitude);
        
        obsIdx = true([1, nvars]);
        leftPad = 100;
        rightPad = pertLength;
        [est, ~, ~, ~, ~, tableResults] = ...
            PerturbationBaseExperiment(noisyData, mat, numTrials, preprocfn, ...
                obsIdx, pertIdx, pertTimes, leftPad, rightPad, method, corrThresh);

        parsave(sprintf('%s/dataLog.mat', currExpPath), noisyData, pertLength, pertTimes, mat, K);
        
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
numRerun = sum(reshape(numRerun, [noiseMagnitudesLength, numMats]), 4);

% Save experiment results
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', 'accLog', 'numRerun');


%% Plot Results

% Show number of simulations that were skipped.
figure(1)
plot(noiseMagnitudes, numRerun)
xlabel('Magnitude of Noise')
ylabel('Number of Simulations Rerun')


% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 2);
figure(2)
plot(noiseMagnitudes, aveAccuracies)
xlabel('Magnitude of Noise')
ylabel('Average Accuracy over Simulations')


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 4);
figure(3)
plot(noiseMagnitudes, aveTPR)
xlabel('Magnitude of Noise')
ylabel('Average TPR over Simulations')


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 4);
figure(4)
plot(noiseMagnitudes, aveFPR)
xlabel('Magnitude of Noise')
ylabel('Average FPR over Simulations')
