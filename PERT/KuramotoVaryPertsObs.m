clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')
addpath('./BAgraph_dir/')
addpath('./kmeans_opt/')

warning('off', 'stats:kmeans:MissingDataRemoved')

%% Initialize Parameters

% Number of network nodes
nvars = 20;

% Probability of ER network edge connections
prob = 0.5;

% Gaussian noise function
noiseVar = 0;
noisefn = @(data) WhiteGaussianNoise(data, noiseVar);

% Delta t
deltat = 0.1;

% Initial conditions
pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);
cfn = @(n) constfn(n, damping);

% Connection strength
K = 50;

% Perturbation force for oscillators
pertForce = 50;

% Width of moving variance window
movvarWidth = 20;

% Threshold for mean variance algorithm
meanThresh = 0;

% Padding for window
pad = 0;

% Require condition that all perturbed nodes must be observed
pertIsObs = false;

% Number of experimental trials
numTrials = 10;

method = 'meanvar';

% Make directory to hold results files if one does not already exist
expName = sprintf('EXP(nvars%d_prob%.2f_K%.2f_damping%.2f_pertf%.2f)', nvars, prob, K, damping, pertForce);
expPath = sprintf('KuramotoExperiments/%s', expName);
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


%% Simulate Data with Varying Numbers of Perturbations
dataLog = cell(nvars, nvars, numTrials);
trueMats = cell(nvars, nvars, numTrials);
dataParamsLog = cell(nvars, nvars, numTrials, 4);

fprintf('Simulate Data:\n')
for numObs = nvars:nvars
    numObs
    
    % If we can only perturb the nodes we can observe, then the maximum
    % number of nodes we can perturb is the number of observed nodes.
    maxNumPerts = nvars;
    if pertIsObs
        maxNumPerts = numObs;
    end
    
    for numPerts = maxNumPerts:maxNumPerts
        numPerts
        
        % Randomly choose which nodes we are allowed to observe.
        inds = randsample(nvars, numObs);
        obsIdx = false([1, nvars]);
        obsIdx(inds) = 1;
        
        % If we can only perturb the nodes we can observe, choose the
        % nodes to perturb out of the nodes we can observe.
        if pertIsObs
            if nnz(obsIdx) == 1
                % Randsample doesn't work in the edge case where obsIdx
                % contains one element.
                pertIdx = find(obsIdx);
            else
                pertIdx = randsample(find(obsIdx), numPerts);
            end
        else
            pertIdx = randsample(nvars, numPerts);
        end
        
        trial = 1;
        while trial <= numTrials
            % Build up network connectivity
            mat = MakeNetworkER(nvars, prob, true);

            % Create the timespan for the simulation.
            waitTime = 10;

            endtime = waitTime * (numPerts + 1);
            nobs = round(endtime / deltat);
            tSpan = linspace(0, endtime, nobs);

            % Build up forcing function.
            times = round(linspace(0, nobs, numPerts+2));
            pertTimes = times(2:end-1);
            pertLength = round(nobs/(5*(numPerts+1)));

            forcingFunc = zeros([nvars, length(tSpan)]);
            for k=1:numPerts
                forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
            end

            % Generate data with forced perturbations.
            data = GenerateKuramotoData(mat, tSpan, 1, K, pfn, wfn, cfn, forcingFunc);
            noisyData = noisefn(data);

            dataLog{numObs, numPerts, trial} = noisyData;
            trueMats{numObs, numPerts, trial} = mat;
            dataParamsLog{numObs, numPerts, trial, 1} = obsIdx;
            dataParamsLog{numObs, numPerts, trial, 2} = pertIdx;
            dataParamsLog{numObs, numPerts, trial, 3} = pertTimes;
            dataParamsLog{numObs, numPerts, trial, 4} = pertLength;

            trial = trial + 1;

            % Save experiment simulated data, connectivity matrices, and parameters.
            save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
            save(sprintf('%s/trueMats.mat', expPath), 'trueMats');
            save(sprintf('%s/dataParamsLog.mat', expPath), 'dataParamsLog');
        end
    end
end


%% Evaluate Algorithm on Data for Varying Numbers of Perturbations and Observations
truePertOrdersLog = cell(nvars, nvars, numTrials);
predPertOrdersLog = cell(nvars, nvars, numTrials);
predMats = cell(nvars, nvars, numTrials);
accuracyLog = nan(nvars, nvars, numTrials);

fprintf('Run Algorithm:\n')
for numObs = nvars:-1:nvars
    numObs
    for numPerts = nvars:-1:nvars
        numPerts
        
        trial = 1;
        while trial <= numTrials
            data = dataLog{numObs, numPerts, trial};
            if isempty(data)
                trial = trial + 1;
                continue
            end
            
            mat = trueMats{numObs, numPerts, trial};
            obsIdx = dataParamsLog{numObs, numPerts, trial, 1};
            pertIdx = dataParamsLog{numObs, numPerts, trial, 2};
            pertTimes = dataParamsLog{numObs, numPerts, trial, 3};
            pertLength = dataParamsLog{numObs, numPerts, trial, 4};
            
            % Select only the subset of nodes that we can observe.
            observedData = data(obsIdx, :);
            [Aprob, predPertOrders] = CreateProbabilityMatrix(observedData, ...
                                                            pertIdx, obsIdx, pertTimes, ...
                                                            pertLength, method, meanThresh, ...
                                                            pad, movvarWidth);
            predPertOrdersLog{numObs, numPerts, trial} = predPertOrders;
            
            % Get the true perturbation orders.
            truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
            truePertOrdersLog{numObs, numPerts, trial} = truePertOrders;
            
            % Get the network reconstruction our algorithm produces.
            predMat = Aprob > 0.5;
            predMats{numObs, numPerts, trial} = predMat;
            
            % Compute the accuracy of this reconstruction (fraction of edges
            % correctly predicted) and store it in the scoring matrix.
            acc = 1 - nnz(predMat(obsIdx, :) - mat(obsIdx, :)) / (nvars * numObs);
            accuracyLog(numObs, numPerts, trial) = acc;
            
            % Save experiment results.
            save(sprintf('%s/truePertOrdersLog.mat', expPath), 'truePertOrdersLog');
            save(sprintf('%s/predPertOrdersLog.mat', expPath), 'predPertOrdersLog');
            save(sprintf('%s/predMats.mat', expPath), 'predMats');
            save(sprintf('%s/accuracyLog.mat', expPath), 'accuracyLog');
            
            trial = trial + 1;
        end
    end
end

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = mean(accuracyLog, 3);
clims = [0, 1];
imagesc(aveAccuracies, clims)
colorbar
xlabel('Number of Perturbed Nodes')
ylabel('Number of Observed Nodes')
