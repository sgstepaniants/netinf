clear all; close all; clc;
addpath('../DataScripts/SimulateData/')

expNum = 'Paper1';

% Check that directory with experiment data exists
expName = sprintf('EXP%s', expNum);
expPath = sprintf('../KuramotoExperiments/%s', expName);

if exist(expPath, 'dir') ~= 7
    error('Data not found: %s', expPath)
end

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

load(sprintf('%s/params.mat', expPath))
load(sprintf('%s/dataLog.mat', expPath))
load(sprintf('%s/trueMats.mat', expPath))

% Run Granger Causality to infer network connections.
preprocfn = @(data) cos(data);
save(sprintf('%s/expParams.mat', resultPath), 'preprocfn')

[predMats, tableResults] = GrangerBaseExperiment(dataLog, trueMats, preprocfn);

% Create a confusion matrix for network predictions.
confusionMat = ConfusionMatrix(trueMats, predMats)

save(sprintf('%s/predMats.mat', resultPath), 'predMats');
save(sprintf('%s/tableResults.mat', resultPath), 'tableResults');
save(sprintf('%s/confusionMat.mat', resultPath), 'confusionMat');
