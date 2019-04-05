clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')
addpath('./BAgraph_dir/')
addpath('./kmeans_opt/')

nvars = 5;
prob = 0.5;
spring = 0.1;
damping = 0.3;
pertForce = 30;

expNum = sprintf('EXP(nvars%d_prob%.2f_spring%.2f_damping%.2f_pertf%.2f)', nvars, prob, spring, damping, pertForce);
expPath = sprintf('HarmonicExperiments/%s', expNum);

if exist(expPath, 'dir') ~= 7
    error('Data not found: %s', expPath)
end

resultPath = sprintf('%s/PertResults', expPath);
if exist(resultPath, 'dir') ~= 7
    mkdir(resultPath)
else
    m=input(sprintf('%s\n already exists, would you like to continue and overwrite these results (Y/N): ', resultPath),'s');
    if upper(m) == 'N'
       return
    end
end

% Threshold for correlation algorithm
corrThresh = 0.5;

% Padding for window
pad = 100;

% Method for building probability matrix
method = 'corr';

%% Evaluate Algorithm on Data for Varying Numbers of Perturbations and Observations
load(sprintf('%s/params.mat', expPath))

predMats = nan(nvars, nvars, numMats, numTrials, nvars, nvars);
accuracyLog = nan(nvars, nvars, numMats, numTrials);

fprintf('Run Algorithm:\n')
for numObs = nvars:-1:1
    numObs
    for numPerts = nvars:-1:1
        numPerts
        
        try
            load(sprintf('%s/numobs%d/numperts%d/dataLog.mat', expPath, numObs, numPerts), 'dataLog');
            load(sprintf('%s/numobs%d/numperts%d/trueMats.mat', expPath, numObs, numPerts), 'trueMats');
            load(sprintf('%s/numobs%d/numperts%d/dataObsIdx.mat', expPath, numObs, numPerts), 'dataObsIdx');
            load(sprintf('%s/numobs%d/numperts%d/dataPertIdx.mat', expPath, numObs, numPerts), 'dataPertIdx');
            load(sprintf('%s/numobs%d/numperts%d/dataPertTimes.mat', expPath, numObs, numPerts), 'dataPertTimes');
            load(sprintf('%s/numobs%d/numperts%d/dataPertLength.mat', expPath, numObs, numPerts), 'dataPertLength');
        catch
            continue
        end
        
        matCount = 1;
        while matCount <= numMats
            data = dataLog{matCount};
            if isempty(data)
                matCount = matCount + 1;
                continue
            end
            
            mat = squeeze(trueMats(matCount, :, :));
            obsIdx = dataObsIdx(matCount, :);
            pertIdx = dataPertIdx(matCount, :);
            pertTimes = dataPertTimes(matCount, :);
            pertLength = dataPertLength(matCount);
            
            for trial = 1:numTrials
                % Select only the subset of nodes that we can observe.
                observedData = data(obsIdx, :, :);
                [AprobHist, predPertOrders] = CreateProbabilityMatrix(observedData, ...
                                                                pertIdx, obsIdx, pertTimes, ...
                                                                pertLength, method, corrThresh, pad);

                % Get the true perturbation orders.
                %truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);

                % Get the network reconstruction our algorithm produces.
                AprobLast = AprobHist(:, :, end);
                predMat = double(AprobLast > 0.5);
                predMat(isnan(AprobLast)) = NaN;
                predMats(numObs, numPerts, matCount, trial, :, :) = predMat;

                % Compute the accuracy of this reconstruction (fraction of edges
                % correctly predicted) and store it in the scoring matrix.
                acc = nnz((predMat == mat) .* eye(nvars)) / (numObs^2-numObs);
                accuracyLog(numObs, numPerts, matCount, trial) = acc;
                
                tpr = 
                fpr = nnz(predMat - mat == 1);

                % Save experiment results.
                save(sprintf('%s/PertResults/predMats.mat', expPath), 'predMats');
                save(sprintf('%s/PertResults/accuracyLog.mat', expPath), 'accuracyLog');
            end
            
            matCount = matCount + 1;
        end
    end
end

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = mean(mean(accuracyLog, 4), 3);
clims = [0, 1];
imagesc(aveAccuracies, clims)
colorbar
xlabel('Number of Perturbed Nodes')
ylabel('Number of Observed Nodes')
