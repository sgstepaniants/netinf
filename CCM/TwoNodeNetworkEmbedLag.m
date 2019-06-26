clear all; close all; clc;
run('../mvgc_v1.0/startup.m')
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')

expNum = 'ConfusionMatrix_Size2';
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../HarmonicExperiments/%s', expName);

maxLag = 100;
maxEmb = 10;
fnnPercentThresh = 5;

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


%% Run mdembedd to find optimal time lags and embedding dimension
load(sprintf('%s/params.mat', expPath), '-regexp', '^(?!expNum$|expName$|expPath$|resultPath$|preprocfn$).')

maxTrials = 1;

tauLog = nan(maxTrials, numMats);
ELog = nan(nvars, maxTrials, numMats);
fnnPercentLog = nan([nvars, maxTrials, numMats, maxEmb]);

load(sprintf('%s/dataLog.mat', expPath));
for j = 1 : maxTrials
    j
    for k = 1 : numMats;
        data = dataLog(:, :, j);

        tau = mdDelay(data.', 'maxLag', maxLag, 'plottype', 'none');
        tauLog(j, k) = tau;
        
        for n = 1 : nvars
            fnnPercent = mdFnn(data(n, :).', round(tau), 'maxEmb', maxEmb, 'doPlot', 0);
            fnnPercentLog(n, j, k, :) = fnnPercent;
            
            E = find(fnnPercent < fnnPercentThresh, 1, 'first');
            if isempty(E)
                E = NaN;
            end
            ELog(n, j, k) = E;
        end
    end
end

tauNone = nanmean(nanmean(tauLog(:, 1:100)))
tau1causes2 = nanmean(nanmean(tauLog(:, 101:200)))
tau2causes1 = nanmean(nanmean(tauLog(:, 201:300)))
tauBoth = nanmean(nanmean(tauLog(:, 301:400)))

ENone = nanmean(nanmean(nanmean(ELog(:, :, 1:100))))
E1causes2 = nanmean(nanmean(nanmean(ELog(:, :, 101:200))))
E2causes1 = nanmean(nanmean(nanmean(ELog(:, :, 201:300))))
EBoth = nanmean(nanmean(nanmean(ELog(:, :, 301:400))))

save(sprintf('%s/results.mat', resultPath), 'tauLog', 'ELog', 'fnnPercentLog', ...
    'tauNone', 'tau1causes2', 'tau2causes1', 'tauBoth', ...
    'ENone', 'E1causes2', 'E2causes1', 'EBoth');
