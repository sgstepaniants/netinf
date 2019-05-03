clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'Paper3';

% Network size
nvars = 3;

% Initialize masses, positions, and velocities of oscillators.
mfn = @(n) constfn(n, 1); %randfn(n, 0.7, 1.3);
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) zeros([n, 1]);

% Specify the damping constant.
damping = 0.2;
cfn = @(n) constfn(n, damping);

% Define time sampling.
deltat = 0.1; % space between time points
endtime = 25;
nobs = round(endtime / deltat); % number of time points (observations)
tSpan = linspace(0, endtime, nobs);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn  = @(data) WhiteGaussianNoise(data, measParam);

% Specify forcing function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Specify boundary conditions.
bc = 'fixed';

% Probabilities of network connections to try.
probs = 0.5; %0.1:0.1:1;
numProbs = length(probs);

% Connection strengths to try.
strengths = 1; %1:1:10;
numStrengths = length(strengths);

% Number of matrices to try for each probability and connection strength
% combination.
numMats = 10;

% Number of simulation trials.
numTrials = 100;

% Spectral radius threshold for MVGC toolbox.
rhoThresh = 0.99;


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
resultPath = sprintf('%s/GCResults', expPath);
if exist(resultPath, 'dir') ~= 7
    mkdir(resultPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
end

%% Generate Data and Run Granger Causality Experiments

% Create random connectivity matrices and simulate oscillator trajectories.
dataLog = nan(nvars, nobs, numTrials, numMats, numProbs, numStrengths);
trueMats = nan(nvars, nvars, numMats, numProbs, numStrengths);
Ks = nan(nvars+2, nvars+2, numMats, numProbs, numStrengths);

% Run Granger Causality to infer network connections.
reps = 1;
freq = 1;
preprocfn = @(data) standardize(data);
save(sprintf('%s/expParams.mat', resultPath), 'reps', 'freq', 'preprocfn')

predMats = nan(nvars, nvars, numMats, numProbs, numStrengths);
est = nan(nvars, nvars, reps, numMats, numProbs, numStrengths);
numRerun = zeros(numProbs, numStrengths);
tprLog = nan(numProbs, numStrengths, numMats);
fprLog = nan(numProbs, numStrengths, numMats);
accuracyLog = nan(numProbs, numStrengths, numMats);
tableResultsLog = repmat(struct([]), [numProbs, numStrengths, numMats]);

for j = 1 : numProbs
    prob = probs(j)
    
    for k = 1 : numStrengths
        strength = strengths(k)
        
        l = 1;
        while l <= numMats
            % Create adjacency matrices.
            mat = MakeNetworkER(nvars, prob, true);
            K = MakeNetworkTriDiag(nvars+2, false);
            K(2:nvars+1, 2:nvars+1) = mat;
            K = strength .* K;
            
            % If any nodes in the network are not connected to the walls or
            % the eigenvalues of the system have positive real parts, don't
            % use this network.
            [~, amplitudes] = checkHarmonicMat(K, damping);
            if any(amplitudes > 0)
                numRerun(j, k) = numRerun(j, k) + 1;
                continue
            end
            
            % Simulate oscillator trajectories.
            data = GenerateHarmonicData(nvars, tSpan, ...
                    numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
            noisyData = noisefn(data);
            
            dataObsIdx = true([1, nvars]); % default parameter
            [predMat, estm, tableResults] = GrangerBaseExperiment(noisyData, ...
                    mat, numTrials, reps, preprocfn, freq, '', dataObsIdx, rhoThresh);
            if isnan(predMat)
                numRerun(j, k) = numRerun(j, k) + 1;
                continue
            end
            
            dataLog(:, :, :, l, j, k) = noisyData;
            trueMats(:, :, l, j, k) = mat;
            Ks(:, :, l, j, k) = K;
            
            predMats(:, :, l, j, k) = predMat;
            est(:, :, :, l, j, k) = estm;
            
            tableResultsLog(j, k, l).tpr = tableResults.tpr;
            tableResultsLog(j, k, l).tprVoting = tableResults.tprVoting;
            tableResultsLog(j, k, l).fpr = tableResults.fpr;
            tableResultsLog(j, k, l).fprVoting = tableResults.fprVoting;
            tableResultsLog(j, k, l).acc = tableResults.acc;
            tableResultsLog(j, k, l).accVoting = tableResults.accVoting;
            tableResultsLog(j, k, l).diagnostics = tableResults.diagnostics;
            
            tprLog(j, k, l) = tableResults.tprVoting;
            fprLog(j, k, l) = tableResults.fprVoting;
            accuracyLog(j, k, l) = tableResults.accVoting;
            
            l = l + 1;
        end
    end
end

% Save experiment simulated data and connectivity matrices.
save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
save(sprintf('%s/trueMats.mat', expPath), 'trueMats');
save(sprintf('%s/Ks.mat', expPath), 'Ks');

% Save experiment results
save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/est.mat', resultPath), 'est');
save(sprintf('%s/tableResultsLog.mat', resultPath), 'tableResultsLog');
save(sprintf('%s/tprLog.mat', resultPath), 'tprLog');
save(sprintf('%s/fprLog.mat', resultPath), 'fprLog');
save(sprintf('%s/accuracyLog.mat', resultPath), 'accuracyLog');


%% Plot Results

% Show number of simulations that were skipped.
figure(1)
imagesc(numRerun)
colorbar
title('Number of Simulations Rerun by Our Analysis')
xlabel('Connection Strength')
ylabel('Connection Probability')

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accuracyLog, 3);
figure(2)
clims = [0, 1];
imagesc(aveAccuracies, clims)
colorbar
title('Average Accuracy over Simulations')
xlabel('Connection Strength')
ylabel('Connection Probability')

% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 3);
figure(3)
clims = [0, 1];
imagesc(aveTPR, clims)
colorbar
title('Average TPR over Simulations')
xlabel('Connection Strength')
ylabel('Connection Probability')

% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 3);
figure(4)
clims = [0, 1];
imagesc(aveFPR, clims)
colorbar
title('Average FPR over Simulations')
xlabel('Connection Strength')
ylabel('Connection Probability')
