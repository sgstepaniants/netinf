clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

nvars = 5;

expNum = sprintf('VaryTimeStrengths_Size%d', nvars);

% Simulation endtimes
endtimes = 5 : 5 : 25;
numEndtimes = length(endtimes);

% Connection strengths
strengths = 3 : 10;
numStrengths = length(strengths);

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

% Specify noise and prepocessing for data.
measParam = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, measParam);

% Delta t
deltat = 0.1;

% Preprocessing function for data.
preprocfn = @(data) cos(data);

% Probabilities of network connections.
prob = 0.5;

% Number of matrices to average results over.
numMats = 100;

% Number of experimental trials
numTrials = 100;

rhoThresh = 1;

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);
if exist(expPath, 'dir') == 7
    %m=input(sprintf('%s\n already exists, would you like to continue and overwrite this data (Y/N): ', expPath),'s');
    %if upper(m) == 'N'
    %    return
    %end
    rmdir(expPath, 's')
end
mkdir(expPath)

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/GCResults', expPath);
if exist(resultPath, 'dir') == 7
    %m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    %if upper(m) == 'N'
    %   return
    %end
    rmdir(resultPath, 's')
end
mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments

% Run PCI to infer network connections.
predMats = cell(numEndtimes, numStrengths, numMats);
tprLog = nan(numEndtimes, numStrengths, numMats);
fprLog = nan(numEndtimes, numStrengths, numMats);
accLog = nan(numEndtimes, numStrengths, numMats);
save(sprintf('%s/results.mat', resultPath), 'predMats', 'tprLog', 'fprLog', 'accLog');

% Number of parallel processes
M = 12;
c = progress(numEndtimes * numStrengths * numMats);
parfor (idx = 1 : numEndtimes * numStrengths * numMats, M)
    [j, k, m] = ind2sub([numEndtimes, numStrengths, numMats], idx);
    fprintf('endtime: %d, strength: %d\n', j, k)

    currExpPath = sprintf('%s/endtime%d/strength%d/mat%d', expPath, j, k, m);
    if exist(sprintf('%s/dataLog.mat', currExpPath), 'file') ~= 2
        mkdir(currExpPath)
    else
        continue
    end

    % Count the number of iterations done by the parfor loop
    c.count();

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
            continue
        end

        parSave.parDataSave(sprintf('%s/dataLog.mat', currExpPath), noisyData, mat);
        presults = load(sprintf('%s/results.mat', resultPath));
        parSave.parResultsSave(sprintf('%s/results.mat', resultPath), j, k, m, results, est,...
            tableResults.tpr, tableResults.fpr, tableResults.acc);
        break
    end
end

results = load(sprintf('%s/results.mat', resultPath));
predMats = results.predMats;
tprLog = results.tprLog;
fprLog = results.fprLog;
accLog = results.accLog;


%% Plot Results

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = nanmean(accLog, 3);
figure(2)
clims = [0, 1];
imagesc(aveAccuracies, clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average Accuracy over Simulations')
xlabel('Connection Strength')
ylabel('Endtime')
set(gca, 'XTick', strengths)
set(gca, 'YTick', endtimes)
%set(gca, 'TickLength', [0 0])


% Show average TPR for each number of perturbations and
% observations.
aveTPR = nanmean(tprLog, 3);
figure(3)
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
figure(4)
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
