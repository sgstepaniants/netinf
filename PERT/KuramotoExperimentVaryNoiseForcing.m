clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts')
addpath('../DataScripts/SimulateData')
addpath('../DataScripts/SimulateData/InitFunctions')

expNum = 'VaryNoiseForcing';

% Network size
nvars = 10;

% Noise magnitudes
noiseMagnitudes = 0 : 0.2 : 2;
noiseMagnitudesLength = length(noiseMagnitudes);

% Forcing magnitudes
forces = 10 : 10 : 100;
numForces = length(forces);

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
predMats = cell(1, noiseMagnitudesLength * numForces * numMats);
tprLog = nan(1, noiseMagnitudesLength * numForces * numMats);
fprLog = nan(1, noiseMagnitudesLength * numForces * numMats);
accLog = nan(1, noiseMagnitudesLength * numForces * numMats);

parsave = @(fname, noisyData, pertIdx, obsIdx, pertLength, pertTimes, mat)...
            save(fname, 'noisyData', 'pertIdx', 'obsIdx', 'pertLength', 'pertTimes', 'mat');

% Number of parallel processes
M = 25;
c = progress(noiseMagnitudesLength * numForces * numMats);
parfor (idx = 1 : noiseMagnitudesLength * numForces * numMats, M)
    [j, k, m] = ind2sub([noiseMagnitudesLength, numForces, numMats], idx);
    fprintf('noise magnitude: %d, force: %d\n', j, k)
    
    currExpPath = sprintf('%s/noise%d/force%d/mat%d', expPath, j, k, m);
    if exist(currExpPath, 'dir') ~= 7
        mkdir(currExpPath)
    end
    
    % Count the number of iterations done by the parfor loop
    c.count();

    noiseMagnitude = noiseMagnitudes(j);
    force = forces(k);

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


predMats = reshape(predMats, noiseMagnitudesLength, numForces, numMats);
tprLog = reshape(tprLog, [noiseMagnitudesLength, numForces, numMats]);
fprLog = reshape(fprLog, [noiseMagnitudesLength, numForces, numMats]);
accLog = reshape(accLog, [noiseMagnitudesLength, numForces, numMats]);

% Save experiment results
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', 'accLog');


%% Plot Results

% Show average accuracies
aveAccuracies = nanmean(accLog, 3);
figure(1)
clims = [0, 1];
imagesc(reshape(aveAccuracies, [noiseMagnitudesLength, numForces]), clims)
set(gca,'YDir','normal')
colormap jet
colorbar
title('Average Accuracy')
xlabel('Force')
ylabel('Noise')
set(gca, 'XTick', forces)
set(gca, 'YTick', noiseMagnitudes)
%set(gca,'TickLength', [0 0])


% Show average TPR
aveTPR = nanmean(tprLog, 3);
figure(2)
imagesc(reshape(aveTPR, [noiseMagnitudesLength, numForces]))
set(gca,'YDir','normal')
colormap jet
colorbar
title('Average TPR')
xlabel('Force')
ylabel('Noise')
set(gca, 'XTick', forces)
set(gca, 'YTick', noiseMagnitudes)
%set(gca,'TickLength', [0 0])


% Show average FPR
aveFPR = nanmean(fprLog, 3);
figure(3)
imagesc(reshape(aveFPR, [noiseMagnitudesLength, numForces]))
set(gca,'YDir','normal')
colormap jet
colorbar
title('Average FPR')
xlabel('Force')
ylabel('Noise')
set(gca, 'XTick', forces)
set(gca, 'YTick', noiseMagnitudes)
%set(gca,'TickLength', [0 0])
