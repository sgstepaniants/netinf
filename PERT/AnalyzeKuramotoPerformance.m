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

% Damping coefficient
damping = 0;

% Initial conditions
pfn = @(n) 2*pi*rand([n, 1]); % uniform [0, 2pi]
wfn = @(n) 2*rand([n, 1]) - ones(n,1); % uniform [-1, 1]
cfn = @(n) constcfn(n, damping);

% Connection strength
K = 50;

% Perturbation force for oscillators
pertForce = 50;

% Width of moving variance window
movvarWidth = 2;

% Threshold for mean variance algorithm
meanThresh = 0;

% Make directory to hold results files if one does not already exist
expName = sprintf('EXP(nvars%d_prob%.2f_K%.2f_damping%.2f_pertf%.2f)', nvars, prob, K, damping, pertForce);
expPath = sprintf('KuramotoExperiments/%s', expName);
if exist(expPath, 'dir') ~= 7
    mkdir(expPath)
end

% Save experiment parameters.
save(sprintf('%s/params.mat', expPath));


%% Siumlate Data with Varying Numbers of Perturbations
numTrials = 10;
dataLog = cell(nvars, numTrials, 5);

fprintf('Simulate Data:\n')
for numPerts = 1:nvars
    numPerts
    
    trial = 1;
    while trial <= numTrials
        % Build up network connectivity
        mat = MakeNetworkER(nvars, prob, true);

        % Create the timespan for the simulation.
        waitTime = 10;

        deltat = 0.1;
        endtime = waitTime * (numPerts + 1);
        nobs = round(endtime / deltat);
        tSpan = linspace(0, endtime, nobs);
        
        % Build up forcing function.
        pertIdx = randsample(nvars, numPerts);
        times = round(linspace(0, length(tSpan), numPerts+2));
        pertTimes = times(2:end-1);
        pertLength = round(nobs/(5*(numPerts+1)));
        
        forcingFunc = zeros([nvars, length(tSpan)]);
        for k=1:numPerts
            forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
        end
        
        % Generate data with forced perturbations.
        data = GenerateKuramotoData(mat, tSpan, 1, K, pfn, wfn, cfn, forcingFunc);
        noisyData = noisefn(data);
        
        dataLog{numPerts, trial, 1} = noisyData;
        dataLog{numPerts, trial, 2} = mat;
        dataLog{numPerts, trial, 3} = pertIdx;
        dataLog{numPerts, trial, 4} = pertTimes;
        dataLog{numPerts, trial, 5} = pertLength;
        
        trial = trial + 1;
        
        % Save experiment simulated data.
        save(sprintf('%s/dataLog.mat', expPath), 'dataLog');
    end
end


%% Evaluate Algorithm on Data for Varying Numbers of Perturbations and Observations
truePredMats = cell(nvars, nvars, numTrials, 2);
accuracyLog = zeros(nvars, nvars, numTrials);

fprintf('Run Algorithm:\n')
for numPerts = nvars:-1:1
    numPerts
    for numObs = nvars:-1:1
        numObs
        for trial = 1:numTrials
            data = dataLog{numPerts, trial, 1};
            mat = dataLog{numPerts, trial, 2};
            pertIdx = dataLog{numPerts, trial, 3};
            pertTimes = dataLog{numPerts, trial, 4};
            pertLength = dataLog{numPerts, trial, 5};
            
            % Randomly choose which nodes we are allowed to observe and
            % hide the unobserved data by setting it to NaN.
            obsIdx = randsample(nvars, numObs);
            
            [Aprob, pertOrders] = MeanVarProbabilityMatrix(data, pertIdx, obsIdx, pertTimes, pertLength, movvarWidth, meanThresh);
            truePertOrders = TruePertOrders(mat, pertIdx, obsIdx);
            
            % Get the network reconstruction our algorithm produces.
            predMat = Aprob > 0.5;
            truePredMats{numPerts, numObs, trial, 1} = mat;
            truePredMats{numPerts, numObs, trial, 2} = predMat;
            
            % Compute the accuracy of this reconstruction (fraction of edges
            % correctly predicted) and store it in the scoring matrix.
            acc = 1 - nnz(predMat(obsIdx, :) - mat(obsIdx, :)) / (nvars * numObs);
            accuracyLog(numObs, numPerts, trial) = acc;
            
            % Save experiment results.
            save(sprintf('%s/results.mat', expPath), 'truePredMats', 'accuracyLog');
        end
    end
end

% Show average accuracies for each number of perturbations and
% observations.
aveAccuracies = mean(accuracyLog, 3);
clims = [0, 1];
imagesc(aveAccuracies, clims)
colorbar
