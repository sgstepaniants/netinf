clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

% Network size
nvars = 10;

% Experiment name
expNum = sprintf('VaryNoise_Size%d', nvars);

% Noise magnitudes
noiseMagnitudes = 2; %0.2 : 0.2 : 2;
noiseMagnitudesLength = length(noiseMagnitudes);

% Initialize masses, positions, and velocities of oscillators.
mfn = @(n) constfn(n, 1);
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

prob = 0.5;

strength = 1;

% Specify forcing function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Specify boundary conditions.
bc = 'fixed';

% Number of matrices to try for each probability and connection strength
% combination.
numMats = 100;

% Number of simulation trials.
numTrials = 10;

% Spectral radius threshold for MVGC toolbox.
rhoThresh = 0.995;


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
resultPath = sprintf('%s/GCResults', expPath);
%if exist(resultPath, 'dir') == 7
%    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
%    if upper(m) == 'N'
%       return
%    end
%    rmdir(resultPath, 's')
%end
%mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

% Run Granger Causality to infer network connections.
preprocfn = @(data) standardize(data);
%save(sprintf('%s/expParams.mat', resultPath), 'preprocfn')

% Run Granger Causality to infer network connections.
predMats = nan(nvars, nvars, noiseMagnitudesLength * numMats);
tprLog = nan(1, noiseMagnitudesLength * numMats);
fprLog = nan(1, noiseMagnitudesLength * numMats);
accLog = nan(1, noiseMagnitudesLength * numMats);
numRerun = zeros(1, noiseMagnitudesLength * numMats);
diagnosticsLog = nan(noiseMagnitudesLength * numMats, 3);

parDataSave = @(fname, noisyData, mat, K)...
            save(fname, 'noisyData', 'mat', 'K');    
parResultsSave = @(fname, est, tpr, fpr, acc, diagnostics)...
            save(fname, 'est', 'tpr', 'fpr', 'acc', 'diagnostics');

% Number of parallel processes
M = 12;
c = progress(noiseMagnitudesLength * numMats);
for (idx = 1 : noiseMagnitudesLength * numMats) %parfor (idx = 1 : noiseMagnitudesLength * numMats, M)
    [j, m] = ind2sub([noiseMagnitudesLength, numMats], idx);
    fprintf('noise: %d\n', j)
    
    c.count();
    
    currExpPath = sprintf('%s/noise%d/mat%d', expPath, j, m);
    if exist(sprintf('%s/dataLog.mat', currExpPath), 'file') ~= 2
        mkdir(currExpPath)
    end
    
    noiseMagnitude = noiseMagnitudes(j);
    
    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);
        K = MakeNetworkTriDiag(nvars+2, false);
        K(2:nvars+1, 2:nvars+1) = mat;
        K = strength .* K;

        % If any eigenvalues of the system have positive real parts, don't
        % use this network.
        [~, amplitudes] = checkHarmonicMat(K, damping);
        if any(amplitudes > 0)
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end

        % Simulate oscillator trajectories.
        data = GenerateHarmonicData(nvars, tSpan, ...
                numTrials, K, pfn, vfn, mfn, cfn, bc, forcingFunc);
        noisyData = WhiteGaussianNoise(data, noiseMagnitude);
        
        dataObsIdx = true([1, nvars]); % default parameter
        [est, tableResults] = GrangerBaseExperiment(noisyData, ...
            mat, preprocfn, dataObsIdx, rhoThresh);
        if isnan(est)
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end

        %parDataSave(sprintf('%s/dataLog.mat', currExpPath), noisyData, mat, K);
        %parResultsSave(sprintf('%s/results.mat', currExpPath), est, ...
        %    tableResults.tpr, tableResults.fpr, tableResults.acc, tableResults.diagnostics);

        predMats(:, :, idx) = est;
        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accLog(idx) = tableResults.acc;
        diagnosticsLog(idx, :) = tableResults.diagnostics;
        break
    end
end

% Reshape data structures
predMats = reshape(predMats, [nvars, nvars, noiseMagnitudesLength, numMats]);
tprLog = reshape(tprLog, [noiseMagnitudesLength, numMats]);
fprLog = reshape(fprLog, [noiseMagnitudesLength, numMats]);
accLog = reshape(accLog, [noiseMagnitudesLength, numMats]);
numRerun = sum(reshape(numRerun, [noiseMagnitudesLength, numMats]), 2);
diagnosticsLog = reshape(diagnosticsLog, [noiseMagnitudesLength, numMats, 3]);

%save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', ...
%    'accLog', 'diagnosticsLog', 'numRerun');


%% Plot Results

% Show number of simulations that were skipped.
figure(1)
plot(noiseMagnitudes, numRerun)
xlabel('Noise')
ylabel('Simulations Rerun')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])


% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 2);
figure(2)
plot(noiseMagnitudes, aveAccuracies)
xlabel('Noise')
ylabel('Accuracy')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 2);
figure(3)
plot(noiseMagnitudes, aveTPR)
xlabel('Noise')
ylabel('TPR')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 2);
figure(4)
plot(noiseMagnitudes, aveFPR)
xlabel('Noise')
ylabel('FPR')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
