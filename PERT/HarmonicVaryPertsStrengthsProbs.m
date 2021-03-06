clear all; close all; clc;
addpath('../DataScripts')
addpath('../DataScripts/SimulateData')
addpath('../DataScripts/SimulateData/InitFunctions')

expNum = 'VaryPertsStrengthsProbs';

% Network size
nvars = 10;

% Number of perturbations
numPertsList = nvars; %1 : nvars;
numPertsLength = length(numPertsList);

% Connection strengths
strengths = 0.5; %0.1 : 0.1 : 1;
numStrengths = length(strengths);

% Connection probabilities
probs = 0.5; %0.1 : 0.1 : 1;
numProbs = length(probs);

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

% Magnitude of forcing in perturbations.
force = 50;

% Number of matrices to average results over.
numMats = 10;

% Number of experimental trials
numTrials = 1;

method = 'corr';

% Threshold for correlation algorithm
corrThresh = 0.2;

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../HarmonicExperiments/%s', expName);
%if exist(expPath, 'dir') == 7
%    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
%    if upper(m) == 'N'
%        return
%    end
%    rmdir(expPath, 's')
%end
%mkdir(expPath)

% Save experiment parameters.
%save(sprintf('%s/params.mat', expPath));

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/PertResults', expPath);
%if exist(resultPath, 'dir') == 7
%    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
%    if upper(m) == 'N'
%       return
%    end
%    rmdir(resultPath, 's')
%end
%mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

% Run PCI to infer network connections.
predMats = cell(1, numPertsLength * numProbs * numStrengths * numMats);
tprLog = nan(1, numPertsLength * numProbs * numStrengths * numMats);
fprLog = nan(1, numPertsLength * numProbs * numStrengths * numMats);
accLog = nan(1, numPertsLength * numProbs * numStrengths * numMats);
numRerun = zeros(1, numPertsLength * numProbs * numStrengths * numMats);

parDataSave = @(fname, noisyData, pertIdx, pertLength, pertTimes, mat, K)...
            save(fname, 'noisyData', 'pertIdx', 'pertLength', 'pertTimes', 'mat', 'K');
parResultsSave = @(fname, est, tpr, fpr, acc, diagnostics)...
            save(fname, 'est', 'tpr', 'fpr', 'acc');

% Number of parallel processes
M = 12;
c = progress(numPertsLength * numProbs * numStrengths * numMats);
for idx = 1 : numPertsLength * numProbs * numStrengths * numMats %parfor (idx = 1 : numPertsLength * numProbs * numStrengths * numMats, M)
    [j, k, l, m] = ind2sub([numPertsLength, numProbs, numStrengths, numMats], idx);
    fprintf('perts: %d, probs: %d, strengths: %d\n', j, k, l)
    
    % Count the number of iterations done by the parfor loop
    c.count();
    
    numPerts = numPertsList(j);
    prob = probs(k);
    strength = strengths(l);
    
    pertIdx = randsample(1 : nvars, numPerts);
    
    currExpPath = sprintf('%s/numPerts%d/prob%d/strength%d/mat%d', expPath, j, k, l, m);
    if exist(sprintf('%s/dataLog.mat', currExpPath), 'file') ~= 2
        %mkdir(currExpPath)
    end
    
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
        obsIdx = true(1, nvars);
        [est, ~, ~, ~, predPertOrders, tableResults] = ...
            PerturbationBaseExperiment(noisyData, mat, numTrials, preprocfn, ...
                obsIdx, pertIdx, pertTimes, leftPad, rightPad, method, corrThresh);
        
        %parDataSave(sprintf('%s/dataLog.mat', currExpPath), noisyData, pertIdx, pertLength, pertTimes, mat, K);
        %parResultsSave(sprintf('%s/results.mat', currExpPath), est, ...
        %     tableResults.tpr, tableResults.fpr, tableResults.acc);
        
        predMats{idx} = est;
        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accLog(idx) = tableResults.acc;
        break
    end
end

predMats = reshape(predMats, numPertsLength, numProbs, numStrengths, numMats);
tprLog = reshape(tprLog, [numPertsLength, numProbs, numStrengths, numMats]);
fprLog = reshape(fprLog, [numPertsLength, numProbs, numStrengths, numMats]);
accLog = reshape(accLog, [numPertsLength, numProbs, numStrengths, numMats]);
numRerun = sum(reshape(numRerun, [numPertsLength, numProbs, numStrengths, numMats]), 4);

% Save experiment results
%save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', 'accLog', 'numRerun');


%% Plot Results

pertInd = 10;

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 4);
figure(2)
clims = [0, 1];
imagesc(reshape(aveAccuracies(pertInd, :, :), [numProbs, numStrengths]), clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
%colorbar
%title('Average Accuracy over Simulations')
%xlabel('Number of Observations')
%ylabel('Number of Perturbations')
%set(gca, 'XTick', numObsList)
%set(gca, 'YTick', numPertsList)


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 4);
figure(3)
clims = [0, 1];
imagesc(reshape(aveTPR(pertInd, :, :), [numProbs, numStrengths]), clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
%colorbar
%title('Average TPR over Simulations')
%xlabel('Number of Observations')
%ylabel('Number of Perturbations')
%set(gca, 'XTick', numObsList)
%set(gca, 'YTick', numPertsList)


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 4);
figure(4)
clims = [0, 1];
imagesc(reshape(aveFPR(pertInd, :, :), [numProbs, numStrengths]), clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
%colorbar
%title('Average FPR over Simulations')
%xlabel('Number of Observations')
%ylabel('Number of Perturbations')
%set(gca, 'XTick', numObsList)
%set(gca, 'YTick', numPertsList)


probInd = 8;
strengthInds = 1:10;
plot(1 : 10, squeeze(aveAccuracies(:, probInd, strengthInds)), 'LineWidth', 4)
ylim([0, 1])
set(gca, 'XTick', [])
set(gca, 'YTick', [])

coeffs = zeros(2, 10, 10);
for k = 1 : 10
    for j = 1 : 10
        coeffs(:, k, j) = polyfit(1:10, squeeze(aveAccuracies(:, k, j)).', 1);
    end
end

surf(squeeze(coeffs(1, :, :)))
