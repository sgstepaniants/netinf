clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts')
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

nvars = 15;

expNum = sprintf('VaryTimeStrengths_Size%d', nvars);

% Connection strengths
strengths = 3;
numStrengths = length(strengths);

% Simulation endtimes
endtimes = 5; %5 %5 : 5 : 25; %max(3, 25 / strengths); %5 : 5 : 25;
numEndtimes = length(endtimes);

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Specify noise and prepocessing for data
measParam = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, measParam);

% Delta t
deltat = 0.1;

% Preprocessing function for data.
preprocfn = @(data) cos(data);

% Probabilities of network connections.
prob = 0.5;

% Number of matrices to average results over.
numMats = 1;

% Number of experimental trials
numTrials = 100;

rhoThresh = 0.995;

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);
%if exist(expPath, 'dir') == 7
%    m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
%    if upper(m) == 'N'
%        return
%    end
%    rmdir(expPath, 's')
%end
mkdir(expPath)

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/GCResults', expPath);
%if exist(resultPath, 'dir') == 7
%    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
%    if upper(m) == 'N'
%       return
%    end
%    rmdir(resultPath, 's')
%end
mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

% Run PCI to infer network connections.
predMats = cell(1, numEndtimes * numStrengths * numMats);
tprLog = nan(1, numEndtimes * numStrengths * numMats);
fprLog = nan(1, numEndtimes * numStrengths * numMats);
accLog = nan(1, numEndtimes * numStrengths * numMats);
numRerun = zeros(1, numEndtimes * numStrengths * numMats);
diagnosticsLog = nan(numEndtimes * numStrengths * numMats, 3);

parsave = @(fname, noisyData, mat)...
            save(fname, 'noisyData', 'mat');

% Number of parallel processes
M = 12;
c = progress(numEndtimes * numStrengths * numMats);
for idx = 1 : numEndtimes * numStrengths * numMats %parfor (idx = 1 : numEndtimes * numStrengths * numMats, M)
    [j, k, m] = ind2sub([numEndtimes, numStrengths, numMats], idx);
    fprintf('endtime: %d, strength: %d\n', j, k)
    
    % Count the number of iterations done by the parfor loop
    c.count();

    currExpPath = sprintf('%s/endtime%d/strength%d/mat%d', expPath, j, k, m);
    if exist(sprintf('%s/dataLog.mat', currExpPath), 'file') ~= 2
        mkdir(currExpPath)
    end

    endtime = endtimes(j);
    strength = strengths(k);
    
    while true
        % Create adjacency matrices.
        mat = MakeNetworkER(nvars, prob, true);

        nobs = round(endtime / deltat);
        tSpan = linspace(0, endtime, nobs);

        % Forcing function
        forcingFunc = zeros([nvars, nobs]);

        % Generate data with forced perturbations.
        data = GenerateKuramotoData(mat, tSpan, numTrials, strength, pfn, wfn, forcingFunc);
        noisyData = noisefn(data);
        
        dataObsIdx = true([1, nvars]); % default parameter
        [est, tableResults] = GrangerBaseExperiment(noisyData, ...
                mat, preprocfn, dataObsIdx, rhoThresh);
        if isnan(est)
            numRerun(idx) = numRerun(idx) + 1;
            continue
        end

        parsave(sprintf('%s/dataLog.mat', currExpPath), noisyData, mat);
        
        predMats{idx} = est;
        tprLog(idx) = tableResults.tpr;
        fprLog(idx) = tableResults.fpr;
        accLog(idx) = tableResults.acc;
        diagnosticsLog(idx, :) = tableResults.diagnostics;
        break
    end
end

% Reshape data structures
predMats = reshape(predMats, numEndtimes, numStrengths, numMats);
tprLog = reshape(tprLog, [numEndtimes, numStrengths, numMats]);
fprLog = reshape(fprLog, [numEndtimes, numStrengths, numMats]);
accLog = reshape(accLog, [numEndtimes, numStrengths, numMats]);
diagnosticsLog = reshape(diagnosticsLog, [numEndtimes, numStrengths, numMats, 3]);
numRerun = sum(reshape(numRerun, [numEndtimes, numStrengths, numMats]), 4);

% Save experiment results
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', ...
    'accLog', 'diagnosticsLog', 'numRerun');


%% Plot Results

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 3);
figure(1)
clims = [0, 1];
imagesc(aveAccuracies, clims)
set(gca,'YDir','normal')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
colormap jet
colorbar
title('Average Accuracy over Simulations')
xlabel('Connection Strength')
ylabel('Endtime')
set(gca, 'XTick', strengths)
set(gca, 'YTick', endtimes)
bestEndtimes = 4.5 * nvars ./ strengths;
hold on; line(1:numStrengths, bestEndtimes / 5, 'Linewidth', 5, 'Color', 'k')


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 3);
figure(2)
clims = [0, 1];
imagesc(aveTPR, clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average TPR over Simulations')
xlabel('Connection Strength')
ylabel('Endtime')
set(gca, 'XTick', strengths)
set(gca, 'YTick', endtimes)
%set(gca, 'TickLength', [0 0])


% Show average FPR for each number of perturbations and
% observations.
aveFPR = nanmean(fprLog, 3);
figure(3)
clims = [0, 1];
imagesc(aveFPR, clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average FPR over Simulations')
xlabel('Connection Strength')
ylabel('Endtime')
set(gca, 'XTick', strengths)
set(gca, 'YTick', endtimes)
%set(gca, 'TickLength', [0 0])
