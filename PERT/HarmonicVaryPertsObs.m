clear all; close all; clc;
addpath('../SimulateData/')
addpath('../SimulateData/InitFunctions/')
addpath('./BAgraph_dir/')
addpath('./kmeans_opt/')

%% Initialize Parameters

% Number of network nodes
nvars = 20;
% Boundary conditions of oscillator system
bc = 'fixed';

% Probability of ER network edge connections
prob = 0.5;

% Gaussian noise function
noiseVar = 0.1;
noisefn = @(data) WhiteGaussianNoise(data, noiseVar);

% Initial conditions and masses
pfn = @(n) randfn(n, -0.5, 0.5);
vfn = @(n) randfn(n, -1, 1);
mfn = @(n) constfn(n, 1);

% Spring constants in oscillator network
spring = 0.1;
% Damping in oscillator network
damping = 0.3;

% Delta t
deltat = 0.1;

% Perturbation force for oscillators
pertForce = 30;

% Threshold for correlation algorithm
corrThresh = 0.5;

% Padding for window
pad = 100;

% Require condition that all perturbed nodes must be observed
pertIsObs = true;

% Number of experimental trials
numTrials = 10;

method = 'corr';

% Make directory to hold results files if one does not already exist
expName = sprintf('EXP(nvars%d_prob%.2f_spring%.2f_damping%.2f_pertf%.2f)', nvars, prob, spring, damping, pertForce);
expPath = sprintf('HarmonicExperiments/%s', expName);
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


%% Siumlate Data with Varying Numbers of Perturbations
dataLog = cell(nvars, nvars, numTrials);
trueMats = cell(nvars, nvars, numTrials);
dataParamsLog = cell(nvars, nvars, numTrials, 4);

fprintf('Simulate Data:\n')
for numObs = 1:nvars
    numObs
    
    % If we can only perturb the nodes we can observe, then the maximum
    % number of nodes we can perturb is the number of observed nodes.
    maxNumPerts = nvars;
    if pertIsObs
        maxNumPerts = numObs;
    end
    
    for numPerts = 1:maxNumPerts
        numPerts
        
        % Randomly choose which nodes we are allowed to observe.
        inds = randsample(nvars, numObs);
        obsIdx = false([1, nvars]);
        obsIdx(inds) = 1;
        
        % If we can only perturb the nodes we can observe, choose the
        % nodes to perturb out of the nodes we can observe.
        if pertIsObs
            pertIdx = randsample(find(obsIdx), numPerts);
        else
            pertIdx = randsample(nvars, numPerts);
        end

        trial = 1;
        while trial <= numTrials
            % Build up network connectivity
            mat = MakeNetworkER(nvars, prob, true);
            K = MakeNetworkTriDiag(nvars + 2, false);
            K(2:nvars+1, 2:nvars+1) = mat;
            K = spring * K;

            % Check if this adjacency matrix has disconnected oscillators.
            G = digraph(mat.');
            distLeft = distances(G, 1);
            distRight = distances(G, nvars);
            disconnectedNodes = find(~isfinite(distLeft) & ~isfinite(distRight));

            % Check if this adjacency matrix is resonant.
            A = K(2:nvars+1, 2:nvars+1);
            A = A - diag(sum(K(2:nvars+1, :), 2));
            lambdas = eig(A);
            amplitudes = real(-damping + sqrt(damping^2 + 4 * lambdas)) / 2;

            % Create the timespan for the simulation.
            eps = 0.1;
            waitTime = ceil(log(eps / min(sqrt(sum(((pertForce * inv(A)).^2))))) / max(amplitudes));

            % If this adjacency matrix is bad, make a new simulation.
            if ~isempty(disconnectedNodes) || any(amplitudes > 0.00001) || waitTime > 500
                continue
            end

            endtime = waitTime * (numPerts + 1);
            nobs = round(endtime / deltat);
            tSpan = linspace(0, endtime, nobs);

            % Build up forcing function.
            times = round(linspace(0, nobs, numPerts+2));
            pertTimes = times(2:end-1);
            pertLength = round(nobs/(10*(numPerts+1)));

            forcingFunc = zeros([nvars, nobs]);
            for k=1:numPerts
                forcingFunc(pertIdx(k), pertTimes(k):pertTimes(k)+pertLength) = pertForce;
            end

            % Generate data with forced perturbations.
            data = GenerateNNCoupledData(nvars, tSpan, 1, K, pfn, vfn, ...
                mfn, @(n)constcfn(n, damping), bc, forcingFunc);
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
for numObs = nvars:-1:1
    numObs
    for numPerts = nvars:-1:1
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
                                                            pertLength, method, corrThresh, pad);
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
