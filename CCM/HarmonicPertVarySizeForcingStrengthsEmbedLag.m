clear all; close all; clc;
run '../mvgc_v1.0/startup.m'
addpath('../DataScripts/SimulateData/')

expNum = 'PertVarySizeForcingStrengths';

% Preprocessing function for data.
preprocfn = @(data) data;

maxLag = 1000;
maxEmb = 20;
fnnPercentThresh = 20;

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);

% Make directory to hold result files if one does not already exist
resultPath = sprintf('%s/MdembeddResults', expPath);
if exist(resultPath, 'dir') == 7
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
    rmdir(resultPath, 's')
end
mkdir(resultPath)


%% Generate Data and Run Granger Causality Experiments
load(sprintf('%s/params.mat', expPath), '-regexp', '^(?!expNum$|expName$|expPath$|resultPath$|preprocfn$).')

% Run mdembedd to find optimal embedding dimension and time lags.
tauLog = nan(1, numSizes * numForces * numStrengths * numMats);
ELog = nan(1, numSizes * numForces * numStrengths * numMats);
fnnPercentLog = nan([maxEmb, numSizes * numForces * numStrengths * numMats]);

% Number of parallel processes
M = 25;
c = progress(numSizes * numForces * numStrengths * numMats);
for idx = 1 : numSizes * numForces * numStrengths * numMats %parfor (idx = 1 : numSizes * numForces * numStrengths * numMats, M)
    [j, k, l, m] = ind2sub([numSizes, numForces, numStrengths, numMats], idx);
    fprintf('size: %d, force: %d, strength: %d\n', j, k, l)
    
    currExpPath = sprintf('%s/size%d/force%d/strength%d/mat%d', expPath, j, k, l, m);
    if exist(currExpPath, 'dir') ~= 7
        continue
    end
    
    % Count the number of iterations done by the parfor loop
    c.count();
    
    dataLog = load(sprintf('%s/dataLog.mat', currExpPath));
    data = dataLog.noisyData;
    
    tau = mdDelay(data.', 'maxLag', maxLag, 'plottype', 'none');
    fnnPercent = mdFnn(data(1, :).', round(tau), 'maxEmb', maxEmb, 'doPlot', 0);
    E = find(fnnPercent < fnnPercentThresh, 1, 'first');
    
    tauLog(idx) = tau;
    ELog(idx) = E;
    fnnPercentLog(:, idx) = fnnPercent;
end

% Reshape data structures
tauLog = reshape(tauLog, [numSizes, numForces, numStrengths, numMats]);
ELog = reshape(ELog, [numSizes, numForces, numStrengths, numMats]);
fnnPercentLog = reshape(fnnPercentLog, [maxEmb, numSizes, numForces, numStrengths, numMats]);

% Save experiment results
save(sprintf('%s/tauLog.mat', resultPath), 'tauLog');
save(sprintf('%s/ELog.mat', resultPath), 'ELog');
save(sprintf('%s/fnnPercentLog.mat', resultPath), 'fnnPercentLog');


%% Plot Results

forceInd = 1;

% Show average taus..
aveTaus = nanmean(tauLog, 4);
figure(2)
clims = [0, maxLag];
imagesc(reshape(aveTaus(:, forceInd, :), [numSizes, numStrengths]), clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average Tau over Simulations')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTick', strengths)
set(gca, 'YTick', networkSizes)
%set(gca, 'TickLength', [0 0])


% Show average Es.
aveEs = nanmean(ELog, 4);
figure(3)
clims = [0, maxEmb];
imagesc(reshape(aveEs(:, forceInd, :), [numSizes, numStrengths]), clims)
set(gca,'YDir','normal')
%set(gca, 'XTick', [])
%set(gca, 'YTick', [])
colormap jet
colorbar
title('Average Es over Simulations')
xlabel('Connection Strength')
ylabel('Network Size')
set(gca, 'XTick', strengths)
set(gca, 'YTick', networkSizes)
%set(gca, 'TickLength', [0 0])
